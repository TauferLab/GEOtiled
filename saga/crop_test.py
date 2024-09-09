import geotiledsaga as gts
import time
import os

TILE_SPLIT_SQRT_MIN = 2
TILE_SPLIT_SQRT_MAX = 10

data_path = '/media/volume/geotiled-saga/crop_test'
gts.set_working_directory(data_path)

results = os.path.join(os.getcwd(), 'crop_results.csv')
f = open(results, 'w')
f.write('tile_count,execution_time\n')
f.close()

# Create a text file with download URLs from a shapefile
gts.fetch_dem(shapefile='CA', txt_file='urls', dataset='30m')

# Download files from the created text file
gts.download_files(download_list='urls', download_folder='dem_tiles')

# Build mosaic from DEMs
gts.build_mosaic(input_folder='dem_tiles', output_file='mosaic', description='Elevation')

# Reproject the mosaic to Projected Coordinate System (PCS) EPSG:26918 - NAD83 UTM Zone 18N
gts.reproject(input_file='mosaic', output_file='elevation', projection='EPSG:26918', cleanup=False)

for i in range(TILE_SPLIT_SQRT_MIN,TILE_SPLIT_SQRT_MAX+1):
    # Compute tile count
    tile_count = i**2

    # Record time of cropping
    start_time = time.time()
    gts.crop_into_tiles(input_file='elevation', output_folder='elevation_tiles', num_tiles=tile_count)
    end_time = time.time()

    # Write results
    formatted_results = str(tile_count) + ',' + str(end_time-start_time) + '\n'
    f = open(results, 'a')
    f.write(formatted_results)
    f.close()