from sklearn.metrics import root_mean_squared_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from osgeo import gdal

import geotiledsaga as gts
import numpy as np

import shutil
import time
import glob
import sys
import os

#################
### FUNCTIONS ###
#################

def generate_elevation_data(working_directory, dataset, roi):
    # Set working directory
    gts.set_working_directory(working_directory)

    # Create a text file with download URLs from a shapefile
    gts.fetch_dem(shapefile=roi, txt_file='urls', dataset=dataset)
    
    # Download files from the created text file
    gts.download_files(download_list='urls', download_folder='dem_tiles')

    # Build mosaic from DEMs
    gts.build_mosaic(input_folder='dem_tiles', output_file='mosaic', description='Elevation')
    
    # Reproject the mosaic to Projected Coordinate System (PCS) EPSG:26918 - NAD83 UTM Zone 18N
    gts.reproject(input_file='mosaic', output_file='elevation', projection='EPSG:26918', cleanup=False)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_test(elevation_tif, working_directory, results_file):
    # Set working directory
    gts.set_working_directory(working_directory)
    
    # Convert GeoTIFF to SDAT
    start_time = time.time()
    gts.__convert_file_format(elevation_tif, os.path.join(os.getcwd(), 'elevation.sdat'), "SAGA")

    # Run the SAGA command
    saga_cmd = ["saga_cmd", "-c=1", "ta_morphometry", "0", "-ELEVATION", os.path.join(os.getcwd(), 'elevation.sgrd'), "-SLOPE", os.path.join(os.getcwd(), 'slope.sgrd')]
    gts.__bash(saga_cmd)

    # Convert the SDAT to GeoTIFF
    gts.__convert_file_format(os.path.join(os.getcwd(), 'slope.sdat'), os.path.join(os.getcwd(), 'slope.tif'), "GTiff")

    # Clock execution time
    total_time = time.time() - start_time

    # Write results to file
    formatted_result = str(total_time) + ',' + str(0) + '\n'
    f = open(results_file, 'a')
    f.write(formatted_result)
    f.close()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def geotiled_nodata_test(tile_count, elevation_tif, results_file):
    # Set working directory
    nodata_path = os.path.join(os.getcwd(), 'nodata_test')
    gts.set_working_directory(nodata_path)

    # Create copy of elevation file within new directory
    shutil.copyfile(elevation_tif, os.path.join(os.getcwd(), 'elevation.tif'))
    
    # Get some data on the split elevation tiles
    elevation_tiles_path = os.path.join(os.getcwd(), 'elevation_tiles')
    elevation_files = sorted(glob.glob(elevation_tiles_path + "/*.tif"))
    for elev_file in elevation_files:
        # Get size of each file
        file_size_mb = os.path.getsize(elev_file) / 1048576
        
        # Get nodata percentage of each file
        matrix_data = gdal.Open(elev_file)
        array_data = matrix_data.ReadAsArray()
        nodata_percentage = (len(array_data[array_data == -999999.0]) / (len(array_data) * len(array_data[0]))) * 100
        with open(results_file, 'a') as sys.stdout:
            print('%s | File Size: %.2f MB | Nodata Percentage: %.2f' % (os.path.basename(elev_file), file_size_mb, nodata_percentage))

    with open(results_file, 'a') as sys.stdout:
        print()
        # Run GEOtiled to compute all terrain parameters for each cropped tile
        gts.compute_geotiled(input_folder='elevation_tiles', parameter_list=['SLP'], num_processes=1, num_cores=1, output_folder_prefix='test')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def geotiled_test(elevation_tif, sequential_tif, storage_directory, results_file, tile_count=4, thread_count=1, core_count=1, ):
    # Set working directory
    gts.set_working_directory(storage_directory)

    # Create copy of elevation file within new directory
    shutil.copyfile(elevation_tif, os.path.join(os.getcwd(), 'elevation.tif'))
    
    # Crop reprojected mosaic into specified number of intermediary tiles
    start_time = time.time()
    gts.crop_into_tiles(input_file='elevation', output_folder='elevation_tiles', num_tiles=tile_count)
    crop_time = time.time() - start_time
    
    # Run GEOtiled to compute all terrain parameters for each cropped tile
    start_time = time.time()
    gts.compute_geotiled(input_folder='elevation_tiles', parameter_list=['SLP'], num_processes=thread_count, num_cores=core_count, output_folder_prefix='test')
    compute_time = time.time() - start_time

    # Build mosaics for each of the computed parameters
    start_time = time.time()
    gts.build_mosaic_buffer(input_folder='test_slope_tiles', output_file='slope')
    mosaic_time = time.time() - start_time

    # Get accuracy results against sequentially generated data
    geotiled_data = gdal.Open(os.path.join(os.getcwd(), 'slope.tif'))
    gt_array_data = geotiled_data.ReadAsArray()
    sequential_data = gdal.Open(sequential_tif)
    s_array_data = sequential_data.ReadAsArray()
    
    rmse = root_mean_squared_error(s_array_data, gt_array_data)
    mse = mean_squared_error(s_array_data, gt_array_data)
    mae = mean_absolute_error(s_array_data, gt_array_data)

    # Format all results to a comma-seperated list
    results = [tile_count,rmse,mse,mae,crop_time,compute_time,mosaic_time,0,0,0]
    formatted_results = ''
    for result in results:
        formatted_results = formatted_results + str(result) + ','
    formatted_results = formatted_results[:len(formatted_results)-1] + '\n' # Remove final comma and add new line

    # Save results to CSV
    f = open(results_file, 'a')
    f.write(formatted_results)
    f.close()

############
### MAIN ###
############

### CONSTANTS ###

TILE_COUNTS = [4,9,16,25,36,49,64,81,100]
THREAD_COUNTS = [2,3,4,6,9,12,18]
CORE_COUNTS = [2,4,8,16,32]

### CONFIG ###

# Get commnand line argumments
data_storage_path = sys.argv[1]

# Set working directory
gts.set_working_directory(data_storage_path)

# Create CSVs for storing test results
sequential_results_csv = os.path.join(os.getcwd(), 'sequential_results.csv')
f = open(sequential_results_csv, 'w')
f.write('compute_time,mem_peak\n')
f.close()

tiling_results_csv = os.path.join(os.getcwd(), 'tiling_results.csv')
f = open(tiling_results_csv, 'w')
f.write('tile_count,rmse,mse,mae,crop_time,compute_time,mosaic_time,crop_mem_peak,compute_mem_peak,mosaic_mem_peak\n')
f.close()

threading_results_csv = os.path.join(os.getcwd(), 'threading_results.csv')
f = open(threading_results_csv, 'w')
f.write('thread_count,rmse,mse,mae,crop_time,compute_time,mosaic_time,crop_mem_peak,compute_mem_peak,mosaic_mem_peak\n')
f.close()

core_count_results_csv = os.path.join(os.getcwd(), 'core_count_results.csv')
f = open(core_count_results_csv, 'w')
f.write('cores_per_process,rmse,mse,mae,crop_time,compute_time,mosaic_time,crop_mem_peak,compute_mem_peak,mosaic_mem_peak\n')
f.close()

# Download and create the input elevation data
elevation_data_path = os.path.join(os.getcwd(), 'elevation_data')
generate_elevation_data(elevation_data_path, '30m', 'TN')
elevation_file = os.path.join(elevation_data_path, 'elevation.tif')
gts.set_working_directory(data_storage_path) # Reset working directory

# Get results for sequential execution
sequential_data_path = os.path.join(os.getcwd(), 'sequential_result')
sequential_test(elevation_file, sequential_data_path, sequential_results_csv)
sequential_file = os.path.join(sequential_data_path, 'slope.tif')
gts.set_working_directory(data_storage_path) # Reset working directory

### TESTS ###

# Run tiling tests
tiling_test_path = os.path.join(os.getcwd(), 'tiling_tests')
os.makedirs(tiling_test_path)
for tile_count in TILE_COUNTS:
    tile_directory = os.path.join(tiling_test_path, str(tile_count) + '_tiles')
    geotiled_test(elevation_file, sequential_file, tile_directory, tiling_results_csv, tile_count=tile_count)
gts.set_working_directory(data_storage_path) # Reset working directory

# Run threading tests
threading_test_path = os.path.join(os.getcwd(), 'threading_tests')
os.makedirs(threading_test_path)
for thread_count in THREAD_COUNTS:
    thread_directory = os.path.join(threading_test_path, str(thread_count) + '_threads')
    geotiled_test(elevation_file, sequential_file, thread_directory, threading_results_csv, tile_count=36, thread_count=thread_count)
gts.set_working_directory(data_storage_path) # Reset working directory

# Run core count tests
core_count_test_path = os.path.join(os.getcwd(), 'core_count_tests')
os.makedirs(core_count_test_path)
for core_count in CORE_COUNTS:
    core_count_directory = os.path.join(core_count_test_path, str(core_count) + '_cores')
    geotiled_test(elevation_file, sequential_file, core_count_directory, core_count_results_csv, core_count=core_count)
gts.set_working_directory(data_storage_path) # Reset working directory