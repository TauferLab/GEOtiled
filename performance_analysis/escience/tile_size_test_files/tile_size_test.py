###############
### IMPORTS ###
###############

import tile_size_functions as tsf
import geotiled_esience
import argparse
import time
import os

############
### MAIN ###
############

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("tile_size", type=int)
parser.add_argument("method", type=str)
parser.add_argument("parameter", type=str)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('tile_size_test_results.csv'):
    file = open('tile_size_test_results.csv', 'w')
    file.write("tile_size,method,parameter,run,compute_time\n")
    file.close()

# Initialize variables
compute_time = 0

# Decide which run to do
if (args.method == 'GDAL') or (args.method == 'SAGA'):
    start_time = time.time()
    tsf.sequential_compute(f"ts{args.tile_size}.tif", args.method, args.parameter)
    compute_time = time.time() - start_time
else:
    start_time = time.time()
    converted_method = 'SAGA' if args.method == 'GEOtiled-SG' else 'GDAL'
    geotiled_esience.crop_and_compute(f"ts{args.tile_size}.tif", args.tile_size, args.tile_size, [args.parameter], compute_method=converted_method)
    geotiled_esience.build_mosaic(input_folder=f"unbuffered_{args.parameter}_tiles", output_file=f"{args.parameter}.tif")
    compute_time = time.time() - start_time

# Write results
file = open('tile_size_test_results.csv', 'a')
file.write(f"\"({args.tile_size},{args.tile_size})\",{args.method},{args.parameter},{args.run},{compute_time}\n")
file.close()