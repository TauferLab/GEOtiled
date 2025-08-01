import tile_size_functions as tsf
import argparse
import geotiled
import time
import os

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("method", type=str)
parser.add_argument("input_file", type=str)
parser.add_argument("tile_size", type=int)
parser.add_argument("parameter", type=str)
parser.add_argument("processes", type=int)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('tile_size_test_results.csv'):
    file = open('tile_size_test_results.csv', 'w')
    file.write("method,tile_size,parameter,run,compute_time\n")
    file.close()

# Decide which run to do
compute_time = 0
if (args.method == 'GDAL') or (args.method == 'SAGA'):
    start_time = time.time()
    tsf.sequential_compute(args.input_file, args.method, args.parameter)
    compute_time = time.time() - start_time
else:
    start_time = time.time()
    converted_method = 'SAGA' if args.method == 'GEOtiled-SG' else 'GDAL'
    geotiled.crop_and_compute(args.input_file, [args.parameter], tile_dimensions=[args.tile_size,args.tile_size], compute_method=converted_method, num_processes=args.processes)
    geotiled.mosaic_rasters(input_folder=f"unbuffered_{args.parameter}_tiles", output_file=f"{args.parameter}.tif")
    compute_time = time.time() - start_time

# Write results
file = open('tile_size_test_results.csv', 'a')
file.write(f"{args.method},\"({args.tile_size},{args.tile_size})\",{args.parameter},{args.run},{compute_time}\n")
file.close()