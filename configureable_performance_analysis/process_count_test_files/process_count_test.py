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
if not os.path.exists('process_count_test_results.csv'):
    file = open('process_count_test_results.csv', 'w')
    file.write("method,tile_size,parameter,processes,run,compute_time\n")
    file.close()

# Compute terrain parameter
start_time = time.time()
converted_method = 'SAGA' if args.method == 'GEOtiled-SG' else 'GDAL'
geotiled.crop_and_compute(args.input_file, args.tile_size, args.tile_size, [args.parameter], compute_method=converted_method, num_processes=args.processes)
geotiled.build_mosaic(input_folder=f"unbuffered_{args.parameter}_tiles", output_file=f"{args.parameter}.tif")
compute_time = time.time() - start_time

# Write results
file = open('process_count_test_results.csv', 'a')
file.write(f"{args.method},\"({args.tile_size},{args.tile_size})\",{args.parameter},{args.processes},{args.run},{compute_time}\n")
file.close()