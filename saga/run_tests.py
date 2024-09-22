import geotiledtests as gt
import pandas as pd
import sys

csv_file = sys.argv[1]

test_parameters = pd.read_csv(csv_file)
df = pd.DataFrame(test_parameters)
df = df.reset_index() # Make sure indexes pair with number of rows

for index, row in df.iterrows():
    storage_path = row['storage_path']
    roi = row['region_of_interest']
    dataset = row['dataset']
    tile_counts = [int(i) for i in row['tile_counts'].split()]
    parameter_list = row['parameter_list'].split()
    process_counts = [int(i) for i in row['process_counts'].split()]
    
    gt.run_full_test(storage_path, roi, dataset, tile_counts, parameter_list, process_counts)