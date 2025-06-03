###############
### IMPORTS ###
###############

import optimization_functions as of
import argparse
import geotiled
import time
import os

############
### MAIN ###
############

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("tile_size", type=int)
parser.add_argument("method", type=str)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('optimization_test_results.csv'):
    file = open('optimization_test_results.csv', 'w')
    file.write(f"tile_size,method,run,crop_time,compute_time,mosaic_time\n")
    file.close()

# Initialize variables
crop_time, compute_time, mosaic_time = 0, 0, 0

# Decide which run to do
if (args.method == 'unoptimized'):
    start_time = time.time()
    of.sequential_crop(f"ts{args.tile_size}.tif", "elevation_tiles", args.tile_size, args.tile_size)
    crop_time = time.time() - start_time

    start_time = time.time()
    of.parallel_compute("elevation_tiles", ['slope'], num_processes=64)
    compute_time = time.time() - start_time

    start_time = time.time()
    of.mosaic_average("slope_tiles", "slope.tif")
    mosaic_time = time.time() - start_time
else:
    start_time = time.time()
    of.parallel_crop(f"ts{args.tile_size}.tif", "elevation_tiles", args.tile_size, args.tile_size, num_processes=64)
    crop_time = time.time() - start_time

    start_time = time.time()
    of.parallel_compute("elevation_tiles", ['slope'], num_processes=64)
    compute_time = time.time() - start_time

    start_time = time.time()
    of.mosaic_concat("slope_tiles", "slope.tif", num_processes=64)
    mosaic_time = time.time() - start_time

# Write results
file = open('optimization_test_results.csv', 'a')
file.write(f"\"({args.tile_size},{args.tile_size})\",{args.method},{args.run},{crop_time},{compute_time},{mosaic_time}\n")
file.close()