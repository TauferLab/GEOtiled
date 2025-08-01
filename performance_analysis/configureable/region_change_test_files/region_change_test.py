import argparse
import geotiled
import time
import os

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("method", type=str)
parser.add_argument("region", type=str)
parser.add_argument("tile_size", type=int)
parser.add_argument("parameter", type=str)
parser.add_argument("processes", type=int)
parser.add_argument("run", type=int)
args = parser.parse_args()

# Create results file
if not os.path.exists('region_change_test_results.csv'):
    file = open('region_change_test_results.csv', 'w')
    file.write("region,method,tile_size,parameter,run,compute_time\n")
    file.close()

compute_time = 0
converted_method = 'SAGA' if args.method == 'GEOtiled-SG' else 'GDAL'
if ((converted_method == 'GDAL') and (args.parameter in ['slope','aspect','hillshade'])) or (converted_method == 'SAGA'):
    start_time = time.time()
    
    # Compute terrain parameters
    geotiled.crop_and_compute(f"{args.region}.tif", [args.parameter], tile_dimensions=[args.tile_size,args.tile_size], compute_method=converted_method, num_processes=args.processes)
    
    # Mosaic files
    if (args.parameter == 'channel_network') or (args.parameter == 'drainage_basins'):
        geotiled.merge_shapefiles(input_folder=f"{args.parameter}_tiles", output_file=f"{args.parameter}.shp")
    else:
        geotiled.mosaic_rasters(input_folder=f"unbuffered_{args.parameter}_tiles", output_file=f"{args.parameter}.tif")
        
    compute_time = time.time() - start_time

# Write results
file = open('region_change_test_results.csv', 'a')
file.write(f"{args.region},{args.method},\"({args.tile_size},{args.tile_size})\",{args.parameter},{args.run},{compute_time}\n")
file.close()