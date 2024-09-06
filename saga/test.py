import pandas as pd
import os

PARAMETERS = ['slope','aspect','hillshade','plan_curvature','profile_curvature','convergence_index']

TILE_SPLIT_SQRT_MIN = 2
TILE_SPLIT_SQRT_MAX = 16

data_path = '/media/volume/geotiled-saga/tile_compute_test_ca'
original_results = '/media/volume/geotiled-saga/tile_compute_test_ca/results.csv'
new_results = '/media/volume/geotiled-saga/tile_compute_test_ca/new_results.csv'

f = open(new_results, 'w')
f.write('tile_count,parameter,avg_comp_time_per_tile,avg_comp_time_std,avg_mem_usage_per_tile,avg_mem_usage_std\n')
f.close()

orig_data = pd.read_csv(original_results)
orig_df = pd.DataFrame(orig_data)

for i in range(TILE_SPLIT_SQRT_MIN,TILE_SPLIT_SQRT_MAX+1):
    tile_count = i**2

    csv_file = os.path.join(data_path, str(tile_count) + '_tiles/individual_results.csv')
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    for param in PARAMETERS:
        df_param = df[df['parameter'] == param]
        orig_df_param = orig_df[orig_df['parameter'] == param]
        orig_df_tile = orig_df_param[orig_df_param['tile_count'] == tile_count]

        formatted_results = str(tile_count) + ',' + param + ',' +  orig_df_tile['average_compute_time_per_tile'].to_string(index=False) + ',' + str(df_param['execution_time'].std()) + ',' + orig_df_tile['average_compute_mem_usage_per_tile'].to_string(index=False) + ',' + str(df_param['peak_mem_usage'].std()) + '\n'
        f = open(new_results, 'a')
        f.write(formatted_results)
        f.close()