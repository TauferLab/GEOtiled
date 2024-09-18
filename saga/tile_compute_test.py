# TODO:
# Add standard deviation computation
# Test sequential execution result
# More dynamic run method (specify tile split count and parameters to compute)

from osgeo import gdal

import geotiledsaga as gts
import geotiledsaga2 as gts2
import pandas as pd

import statistics
import shutil
import time
import glob
import sys
import csv
import os

###############################
### PREPROCESSING FUNCTIONS ###
###############################

def create_memory_log_csv(path, file_name='mem_log'):
    log_file = os.path.join(path, file_name+'.csv')
    f = open(log_file, 'w')
    f.write('mem_usage\n')
    f.close()

    return log_file

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def create_metadata_csv(path, file_name='file_metadata'):
    metadata_file = os.path.join(path, file_name+'.csv')
    f = open(metadata_file, 'w')
    f.write('file_name,file_size,nodata_percentage\n')
    f.close()

    return metadata_file

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_results_to_csv(results_file, results):
    formatted_results = ''
    for result in results:
        formatted_results += str(result) + ','
    formatted_results = formatted_results[:-1] + '\n'
    f = open(results_file, 'a')
    f.write(formatted_results)
    f.close()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def start_memory_tracking(tmux_name, log_script, log_file):
    start_cmd = ['tmux', 'new-session', '-d', '-s', tmux_name, 'python', log_script, log_file]
    gts.__bash(start_cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def end_memory_tracking(tmux_name):
    end_cmd = ['tmux', 'kill-session', '-t', tmux_name]
    gts.__bash(end_cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_peak_memory_usage(log_file):
    # Compute peak memory usage
    mem_data = pd.read_csv(log_file)
    df = pd.DataFrame(mem_data)
    mem_usage = df['mem_usage']
    return (max(mem_usage) - min(mem_usage))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_nodata_percentage_of_file(file, nodata_value):
    matrix_data = gdal.Open(file)
    array_data = matrix_data.ReadAsArray()
    return ((len(array_data[array_data == nodata_value]) / (len(array_data) * len(array_data[0]))) * 100)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def generate_elevation_data(dataset, roi):
    # Download elevation files for a specific region of interest within a given dataset
    gts.fetch_dem(shapefile=roi, dataset=dataset, save_to_txt=False, download=True)

    # Build mosaic from DEMs
    gts.build_mosaic(input_folder='dem_tiles', output_file='mosaic', description='Elevation')
    
    # Reproject the mosaic to Projected Coordinate System (PCS) EPSG:26918 - NAD83 UTM Zone 18N
    gts.reproject(input_file='mosaic', output_file='elevation', projection='EPSG:26918', cleanup=False)

#########################
### TESTING FUNCTIONS ###
#########################

def run_sequential_test(elevation_file, parameter_list, results_file):
    # Copy elevation file over
    sequential_path = os.path.join(os.getcwd(), '1_tile')
    gts.set_working_directory(sequential_path)
    new_elevation_file = os.path.join(os.getcwd(), os.path.basename(elevation_file))
    shutil.copyfile(elevation_file, new_elevation_file)

    gts.compute_geotiled(os.path.pathname(new_elevation_file), 
    
    # Convert elevation data to SDAT
    saga_elevation_file = new_elevation_file.replace(".tif", ".sdat")
    gts.__convert_file_format(new_elevation_file, saga_elevation_file, "SAGA")

    # Run computations
    saga_elevation_file_sgrd = saga_elevation_file.replace(".sdat",".sgrd")
    cmd_base = ["saga_cmd", "-c=1"]
    
    if ("slope" in parameter_list):
        cmd_curv = cmd_base + ["ta_morphometry", "0", "-ELEVATION", saga_elevation_file_sgrd, "-SLOPE", saga_file_paths["slope"]]
        gts.__bash(cmd_curv)
    
    if ("aspect" in parameter_list):
        or ("profile_curvature" in param_list) or ("plan_curvature" in param_list):
        # Add command and elevation param
        cmd_curv = cmd_base + ["ta_morphometry", "0", "-ELEVATION", saga_file_paths["elevation"]]

        # Add other requested parameters
        if "slope" in param_list:
            cmd_curv = cmd_curv + ["-SLOPE", saga_file_paths["slope"]]
        if "aspect" in param_list:
            cmd_curv = cmd_curv + ["-ASPECT", saga_file_paths["aspect"]]
        if "profile_curvature" in param_list:
            cmd_curv = cmd_curv + ["-C_PROF", saga_file_paths["profile_curvature"]]
        if "plan_curvature" in param_list:
            cmd_curv = cmd_curv + ["-C_PLAN", saga_file_paths["plan_curvature"]]

        bash(cmd_curv) # Run
    
    if "hillshade" in param_list:
        # Build command and run
        cmd_shade = cmd_base + ["ta_lighting", "0", "-ELEVATION", saga_file_paths["elevation"], "-SHADE", saga_file_paths["hillshade"]]

        bash(cmd_shade) # Run
    
    if "convergence_index" in param_list:
        # Build command
        cmd_ci = cmd_base + ["ta_morphometry", "1", "-ELEVATION", saga_file_paths["elevation"], "-RESULT", saga_file_paths["convergence_index"]]
    

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_crop_test(elevation_file, tile_counts, results_file):
    for tile_count in tile_counts:
        # Create some paths to store results
        tile_path = os.path.join(os.getcwd(), str(tile_count) + '_tiles')
        elevation_tile_path = os.path.join(tile_path, 'elevation_tiles')
        Path(tile_path).mkdir(parents=True, exist_ok=True)
        
        # Start trackers
        mem_log_file = create_memory_log_csv(os.getcwd())
        start_memory_tracking('crop-memory-tracking', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file)
        start_time = time.time()
        
        # Begin cropping
        gts.crop_into_tiles(input_file=elevation_file, output_folder=elevation_tile_path, num_tiles=tile_count)

        # End trackers
        end_time = time.time()
        start_memory_tracking('crop-memory-tracking')

        # Write results of cropping
        write_results_to_csv(results_file, [tile_count, end_time-start_time, compute_peak_memory_usage(mem_log_file)])

        # Get some metadata related to the elevation tiles
        metadata_file = create_metadata_csv(tile_path, 'elevation_metadata')
        elevation_files = sorted(glob.glob(elevation_tile_path))
        for elev_file in elevation_files:
            file_name = os.path.basename(elev_file)
            file_size = os.path.getsize(elev_file)
            nodata_precentage = get_nodata_percentage_of_file(elev_file, -999999.0)

            # Write metadata to file
            write_results_to_csv(metadata_file, [file_name, end_time-file_size, nodata_precentage])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------        

def run_compute_test():

    return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_mosaic_test():

    
    
    return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_full_test(storage_path, dataset, roi, tile_counts, parameter_list):

    run_crop_test()

    run_compute_test()

    run_mosaic_test()
    
    return

############
### MAIN ###
############

PARAMETERS = ['SLP','ASP','HLD','PLC','PFC','CI']
REGIONS_OF_INTEREST = ['DE','TN','CA']
DATASETS = ['30m']
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
f.write('tile_count,parameter,avg_comp_time_per_tile,avg_mem_usage_per_tile\n')
f.close()

file_data_csv = os.path.join(os.getcwd(), 'file_data.csv')
f = open(file_data_csv, 'w')
f.write('array_size,file_size_bytes,nodata_percentage\n')
f.close()

# Download and reproject the elevation data
generate_elevation_data(dataset, region)
elevation_file = os.path.join(os.getcwd(), 'elevation.tif')

# Get metrics for sequential (1 tile) computation
single_tile_path = os.path.join(os.getcwd(), '1_tile')
gts.set_working_directory(single_tile_path)
Path(os.path.join(os.getcwd(), 'elevation_tiles')).mkdir(parents=True, exist_ok=True)
shutil.copyfile(elevation_file, os.path.join(os.getcwd(), 'elevation_tiles/tile_0000.tif'))

individual_results_file = os.path.join(os.getcwd(), 'individual_results.csv')
f = open(individual_results_file, 'w')
f.write('parameter,file_name,execution_time,peak_mem_usage\n')
f.close()

param_codes = gts.__get_codes("parameter")
for parameter in PARAMETERS: 
    # Compute for each parameter
    param = param_codes[parameter]
    gts.__compute_params_saga(os.path.join(os.getcwd(), 'elevation_tiles/tile_0000.tif'), [param, 1, 'ctt', False])
    
# Store results in main result file
data = pd.read_csv(individual_results_file)
df = pd.DataFrame(data)
parameters = df['parameter'].unique() # Get all unique parameters
for param in parameters:
    ex_time_data = df[df['parameter'] == param]['execution_time']
    mem_usage_data = df[df['parameter'] == param]['peak_mem_usage']

    formatted_results = '1,' + param + ',' + str(ex_time_data) + ',' + str(mem_usage_data) + '\n'
    f = open(results_csv, 'a')
    f.write(formatted_results)
    f.close()

gts.set_working_directory(data_storage_path) # Reset working directory

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

    # Create a log file to track execution times and memory usage
    individual_results_file = os.path.join(os.getcwd(), 'individual_results.csv')
    f = open(individual_results_file, 'w')
    f.write('parameter,file_name,execution_time,peak_mem_usage\n')
    f.close()
    
    # Compute for each parameter
    for parameter in PARAMETERS: 
        gts.compute_geotiled(input_folder='elevation_tiles', parameter_list=[parameter], num_processes=1, num_cores=1, output_folder_prefix='ctt')

    # Average together results for each tile and write to file
    data = pd.read_csv(individual_results_file)
    df = pd.DataFrame(data)

    parameters = df['parameter'].unique() # Get all unique parameters

    for param in parameters:
        ex_time_data = df[df['parameter'] == param]['execution_time']
        mem_usage_data = df[df['parameter'] == param]['peak_mem_usage']

        formatted_results = str(tile_count) + ',' + param + ',' + str(statistics.mean(ex_time_data)) + ',' + str(statistics.mean(mem_usage_data)) + '\n'
        f = open(results_csv, 'a')
        f.write(formatted_results)
        f.close()
    
    gts.set_working_directory(data_storage_path) # Reset working directory
