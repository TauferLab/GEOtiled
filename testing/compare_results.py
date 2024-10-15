from sklearn.metrics import mean_squared_error 
from osgeo import gdal
import geotiled as gt
import numpy as np
import time
import os

def compute_accuracy_stats(original_method_file, old_method_file, new_method_file):
    # Open files and read them in as arrays
    orig_data = gdal.Open(original_method_file)
    orig_array = orig_data.ReadAsArray()
    
    old_data = gdal.Open(old_method_file)
    old_array = old_data.ReadAsArray()

    new_data = gdal.Open(new_method_file)
    new_array = new_data.ReadAsArray()

    # Get sums of data (for absolute difference metric)
    orig_sum, old_sum, new_sum = 0, 0, 0
    for i in orig_array:
        orig_sum += sum(i)
    for i in old_array:
        old_sum += sum(i)
    for i in new_array:
        new_sum += sum(i)

    # Output stats
    print('MSE (Original vs Old):', mean_squared_error(orig_array, old_array))
    print('MSE (Original vs New):', mean_squared_error(orig_array, new_array))
    print('Absolute Differnece (Original vs Old):', abs(orig_sum - old_sum))
    print('Absolute Differnece (Original vs New):', abs(orig_sum - new_sum))

def preprocess_data(region, resolution):
    # Download desired data
    gt.fetch_dem(shapefile=region, dataset=resolution, save_to_txt=False, download=True)

    # Mosaic files and reproject
    gt.build_mosaic(input_folder='dem_tiles', output_file='mosaic.tif', description='Elevation')
    gt.reproject(input_file='mosaic.tif', output_file='elevation.tif', projection='EPSG:26918')

def sequential_compute(input_file, output_file):
    # Compute parameter
    seq_compute_time = time.time()
    dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
    gdal.DEMProcessing(output_file, input_file, processing='slope', options=dem_options)
    
    # Set description of data to parameter name
    dataset = gdal.Open(output_file)
    band = dataset.GetRasterBand(1)
    band.SetDescription('slope')
    dataset = None
    print('Original Compute Time:', time.time()-seq_compute_time)

def old_method_compute(input_file, output_file):
    # Crop
    old_crop_time = time.time()
    gt.crop_into_tiles(input_file=input_file, output_folder='tile_elevation_tiles', num_tiles=100)
    print('Old Crop Time:', time.time()-old_crop_time)

    # Compute
    old_compute_time = time.time()
    gt.compute_geotiled(input_folder='tile_elevation_tiles', parameter_list=['slope'], num_processes=10, use_gdal=True, cleanup=True)
    print('Old Compute Time:', time.time()-old_compute_time)

    # Mosaic
    old_mosaic_time = time.time()
    gt.build_mosaic_buffer(input_folder='slope_tiles', output_file=output_file, cleanup=True)
    print('Old Mosaic Time:', time.time()-old_mosaic_time)

def new_method_compute(input_file, output_file):
    # Crop
    new_crop_time = time.time()
    gt.crop_by_size(input_file=input_file, output_folder='size_elevation_tiles', column_length=640, row_length=820)
    print('New Crop Time:', time.time()-new_crop_time)

    # Compute
    new_compute_time = time.time()
    gt.compute_geotiled(input_folder='size_elevation_tiles', parameter_list=['slope'], num_processes=10, use_gdal=True, cleanup=True)
    print('New Compute Time:', time.time()-new_compute_time)

    # Compute
    new_mosaic_time = time.time()
    gt.mosaic_buffered_tiles(input_folder='slope_tiles', output_file=output_file, cleanup=True)
    print('New Mosaic Time:', time.time()-new_mosaic_time)

# Preprocessing
work_dir = '/media/volume/geotiled-saga/compare_test'
gt.set_working_directory(work_dir)

# Set some variables containing file paths
elev_file = os.path.join(work_dir, 'elevation.tif')
orig_file = os.path.join(work_dir, 'seq_slope.tif')
old_method_file = os.path.join(work_dir, 'tile_slope.tif')
new_method_file = os.path.join(work_dir, 'size_slope.tif')

# Download data
preprocess_data('DE','30m')

# Get sequentially computed data
sequential_compute(elev_file, orig_file)

# Compute using old method of cropping and mosaicking
old_method_compute(elev_file, old_method_file)

# Compute using new method of cropping and mosaicking
new_method_compute(elev_file, new_method_file)

# Compare accuracy of final results of old and new method with sequential compute
compute_accuracy_stats(orig_file, old_method_file, new_method_file)