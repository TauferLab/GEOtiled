import optimization_functions as of
import argparse
import time
import os

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("method", type=str)
parser.add_argument("input_file", type=str)
parser.add_argument("tile_size", type=int)
parser.add_argument("library", type=str)
parser.add_argument("parameter", type=str)
parser.add_argument("processes", type=int)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('optimization_test_results.csv'):
    file = open('optimization_test_results.csv', 'w')
    file.write(f"method,tile_size,library,parameter,run,crop_time,compute_time,mosaic_time\n")
    file.close()

# Initialize variables
crop_time, compute_time, mosaic_time = 0, 0, 0

# Decide which run to do
if (args.method == 'unoptimized'):
    start_time = time.time()
    of.sequential_crop(args.input_file, "elevation_tiles", args.tile_size, args.tile_size)
    crop_time = time.time() - start_time

    start_time = time.time()
    of.parallel_compute("elevation_tiles", args.library, args.parameter, num_processes=args.processes)
    compute_time = time.time() - start_time

    start_time = time.time()
    of.mosaic_average(f"{args.parameter}_tiles", f"{args.parameter}.tif")
    mosaic_time = time.time() - start_time
else:
    start_time = time.time()
    of.parallel_crop(args.input_file, "elevation_tiles", args.tile_size, args.tile_size, num_processes=args.processes)
    crop_time = time.time() - start_time

    start_time = time.time()
    of.parallel_compute("elevation_tiles", args.library, args.parameter, num_processes=args.processes)
    compute_time = time.time() - start_time

    start_time = time.time()
    of.mosaic_concat(f"{args.parameter}_tiles", f"{args.parameter}.tif", num_processes=args.processes)
    mosaic_time = time.time() - start_time

# Write results
file = open('optimization_test_results.csv', 'a')
file.write(f"{args.method},\"({args.tile_size},{args.tile_size})\",{args.library},{args.parameter},{args.run},{crop_time},{compute_time},{mosaic_time}\n")
file.close()