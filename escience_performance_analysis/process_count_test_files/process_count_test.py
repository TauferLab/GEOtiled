###############
### IMPORTS ###
###############

import argparse
import geotiled
import time
import os

############
### MAIN ###
############

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("method", type=str)
parser.add_argument("process", type=int)
parser.add_argument("parameter", type=str)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('process_count_test_results.csv'):
    file = open('process_count_test_results.csv', 'w')
    file.write("method,parameter,processes,run,compute_time\n")
    file.close()

start_time = time.time()
converted_method = 'SAGA' if args.method == 'GEOtiled-SG' else 'GDAL'
geotiled.crop_and_compute("ts2500.tif", 2500, 2500, [args.parameter], compute_method=converted_method, num_processes=args.process)
geotiled.build_mosaic(input_folder=f"unbuffered_{args.parameter}_tiles", output_file=f"{args.parameter}.tif")
compute_time = time.time() - start_time

# Write results
file = open('process_count_test_results.csv', 'a')
file.write(f"{args.method},{args.parameter},{args.process},{args.run},{compute_time}\n")
file.close()