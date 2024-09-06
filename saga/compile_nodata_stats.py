from osgeo import gdal
import pandas as pd
import glob
import os

gdal.UseExceptions()

TILE_COUNTS = [4,16,64,256]

data_path = '/media/volume/geotiled-saga/tile_compute_test_ca'

# Create CSV to write results to
results_csv = os.path.join(data_path, 'nodata_stats.csv')
f = open(results_csv, 'w')
f.write('tile_count,file_name,execution_time,peak_mem_usage,nodata_percentage\n')
f.close()

for tile_count in TILE_COUNTS:
    # Create path to elevation tiles and results CSV
    elevation_file_path = os.path.join(data_path, str(tile_count) + '_tiles/elevation_tiles')
    elevation_files = sorted(glob.glob(elevation_file_path + "/*.tif"))
    individual_results_path = os.path.join(data_path, str(tile_count) + '_tiles/individual_results.csv')

    # Open CSV
    data = pd.read_csv(individual_results_path)
    df = pd.DataFrame(data)

    # Want all rows where parameter is convergence index
    ci_results = df[df['parameter'] == 'convergence_index']

    for file in elevation_files:
        file_name = os.path.basename(file)
        ci_elev = ci_results[ci_results['file_name'] == file_name]

        # Compute nodata percentage for the elevation file
        matrix_data = gdal.Open(file)
        array_data = matrix_data.ReadAsArray()
        nodata_percentage = (len(array_data[array_data == -999999.0]) / (len(array_data) * len(array_data[0]))) * 100
        
        # Write formatted results to file
        formatted_results = str(tile_count) + ',' + file_name + ',' + ci_elev['execution_time'].to_string(index=False) + ',' + ci_elev['peak_mem_usage'].to_string(index=False) + ',' + str(nodata_percentage) + '\n'
        f = open(results_csv, 'a')
        f.write(formatted_results)
        f.close()

# Compute and print out some stats
data = pd.read_csv(results_csv)
df = pd.DataFrame(data)

for tile_count in TILE_COUNTS:
    tile_results = df[df['tile_count'] == tile_count]
    nodata_tiles = tile_results[tile_results['nodata_percentage'] == 100]

    print('Tile Count:', tile_count)
    print('Nodata Tile Count:', nodata_tiles.shape[0])
    print('Total Execution Time:', tile_results['execution_time'].sum())
    print('Execution Time of Nodata Tiles:', nodata_tiles['execution_time'].sum())