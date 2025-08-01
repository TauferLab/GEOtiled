###############
### IMPORTS ###
###############

import geotiled_esience
import argparse
import time
import os

############
### MAIN ###
############

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("region", type=str)
parser.add_argument("method", type=str)
parser.add_argument("parameter", type=str)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('region_change_test_results.csv'):
    file = open('region_change_test_results.csv', 'w')
    file.write("region,method,parameter,run,compute_time\n")
    file.close()

compute_time = 0
converted_method = 'SAGA' if args.method == 'GEOtiled-SG' else 'GDAL'
if ((converted_method == 'GDAL') and (args.parameter in ['slope','aspect','hillshade'])) or (converted_method == 'SAGA'):
    start_time = time.time()
    geotiled_esience.crop_and_compute(f"{args.region}.tif", 2500, 2500, [args.parameter], compute_method=converted_method)
    if (args.parameter == 'channel_network') or (args.parameter == 'drainage_basins'):
        geotiled_esience.merge_shapefiles(input_folder=f"{args.parameter}_tiles", output_file=f"{args.parameter}.shp")
    else:
        geotiled_esience.build_mosaic(input_folder=f"unbuffered_{args.parameter}_tiles", output_file=f"{args.parameter}.tif")
    compute_time = time.time() - start_time

# Write results
file = open('region_change_test_results.csv', 'a')
file.write(f"{args.region},{args.method},{args.parameter},{args.run},{compute_time}\n")
file.close()