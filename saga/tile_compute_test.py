from osgeo import gdal

import geotiledsaga as gts

import shutil
import time
import glob
import sys
import csv
import os

#################
### FUNCTIONS ###
#################

def generate_elevation_data(dataset, roi):
    # Create a text file with download URLs from a shapefile
    gts.fetch_dem(shapefile=roi, txt_file='urls', dataset=dataset)
    
    # Download files from the created text file
    gts.download_files(download_list='urls', download_folder='dem_tiles')

    # Build mosaic from DEMs
    gts.build_mosaic(input_folder='dem_tiles', output_file='mosaic', description='Elevation')
    
    # Reproject the mosaic to Projected Coordinate System (PCS) EPSG:26918 - NAD83 UTM Zone 18N
    gts.reproject(input_file='mosaic', output_file='elevation', projection='EPSG:26918', cleanup=False)

############
### MAIN ###
############

PARAMETERS = ['SLP','ASP','HLD','PLC','PFC','CI']
TILE_SPLIT_SQRT_MIN = 2
TILE_SPLIT_SQRT_MAX = 16

# Get commnand line argumments
data_storage_path = sys.argv[1]
dataset = sys.argv[2]
region = sys.argv[3]

# Set working directory
gts.set_working_directory(data_storage_path)

# Create CSVs for storing test results
results_csv = os.path.join(os.getcwd(), 'results.csv')
f = open(results_csv, 'w')
f.write('tile_count,parameter,average_compute_time_per_tile,average_compute_mem_usage_per_tile\n')
f.close()

file_data_csv = os.path.join(os.getcwd(), 'file_data.csv')
f = open(file_data_csv, 'w')
f.write('array_size,file_size_bytes,nodata_percentage\n')
f.close()

generate_elevation_data(dataset, region)
elevation_file = os.path.join(os.getcwd(), 'elevation.tif')

for i in range(TILE_SPLIT_SQRT_MIN,TILE_SPLIT_SQRT_MAX+1):
    # Compute tile count
    tile_count = i**2

    # Tile crop
    cropped_tiles_directory = os.path.join(os.getcwd(), str(tile_count) + '_tiles')
    gts.set_working_directory(cropped_tiles_directory)
    shutil.copyfile(elevation_file, os.path.join(os.getcwd(), 'elevation.tif'))
    gts.crop_into_tiles(input_file='elevation', output_folder='elevation_tiles', num_tiles=tile_count)

    # Save some tile information
    elevation_files = sorted(glob.glob(os.path.join(os.getcwd(), 'elevation_tiles/*.tif')))
    for file in elevation_files:
        file_name = os.path.basename(file)
        file_size = os.path.getsize(file)
        matrix_data = gdal.Open(file)
        array_data = matrix_data.ReadAsArray()
        nodata_percentage = (len(array_data[array_data == -999999.0]) / (len(array_data) * len(array_data[0]))) * 100
        f = open(file_data_csv, 'a')
        f.write(str([len(array_data[0]), len(array_data)]) + ',' + str(file_size) + ',' + str(nodata_percentage) + '\n')
        f.close()
    
    # Compute for each parameter
    for parameter in PARAMETERS: 
        gts.compute_geotiled(input_folder='elevation_tiles', parameter_list=[parameter], num_processes=1, num_cores=1, output_folder_prefix='ctt')

    # Average together results from each tile and write to file
    individual_results_file = os.path.join(os.getcwd(), 'individual_results.csv')
    count = 0
    with open(individual_results_file, mode='r') as file:
        csvFile = csv.reader(file)
        oldParam = ''
        compute_sum = 0
        mem_sum = 0
        for line in csvFile:
            newParam = line[0]
            if count == 0:
                compute_sum += float(line[2])
                mem_sum += float(line[3])
                oldParam = line[0]
                count += 1
            elif newParam == oldParam:
                # If parameter is still the same on next line, add to total compute sum
                compute_sum += float(line[2])
                mem_sum += float(line[3])
            else:
                # Get average compute for a single tile and save results to file
                compute_avg = compute_sum / tile_count
                mem_avg = mem_sum / tile_count
                f = open(results_csv, 'a')
                f.write(str(tile_count) + ',' + oldParam + ',' + str(compute_avg) + ',' + str(mem_avg) + '\n')
                f.close()

                # Reset sum and set oldParam
                compute_sum = float(line[2])
                mem_sum += float(line[3])
                oldParam = line[0]
        # Handle last line
        compute_avg = compute_sum / tile_count
        mem_avg = mem_sum / tile_count
        f = open(results_csv, 'a')
        f.write(str(tile_count) + ',' + oldParam + ',' + str(compute_avg) + ',' + str(mem_avg) + '\n')
        f.close()
    
    gts.set_working_directory(data_storage_path) # Reset working directory