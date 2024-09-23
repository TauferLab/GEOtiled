from matplotlib import pyplot as plt
from pathlib import Path
from osgeo import gdal

import geotiledsaga as gts
import pandas as pd
import numpy as np 

import shutil
import time
import glob
import os

##########################
### TRACKING FUNCTIONS ###
##########################

def write_results_to_csv(results_file, results):
    formatted_results = ''
    for result in results:
        formatted_results += str(result) + ','
    formatted_results = formatted_results[:-1] + '\n'
    f = open(results_file, 'a')
    f.write(formatted_results)
    f.close()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def start_memory_tracking(tmux_name, log_script, log_file):
    start_cmd = ['tmux', 'new-session', '-d', '-s', tmux_name, 'python', log_script, log_file]
    gts.bash(start_cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def end_memory_tracking(tmux_name):
    end_cmd = ['tmux', 'kill-session', '-t', tmux_name]
    gts.bash(end_cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_peak_memory_usage(log_file):
    # Compute peak memory usage
    mem_data = pd.read_csv(log_file)
    df = pd.DataFrame(mem_data)
    mem_usage = df['mem_usage']
    return (max(mem_usage) - min(mem_usage))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_nodata_percentage(file, nodata_value):
    matrix_data = gdal.Open(file)
    array_data = matrix_data.ReadAsArray()
    return ((len(array_data[array_data == nodata_value]) / (len(array_data) * len(array_data[0]))) * 100), (len(array_data) * len(array_data[0]))

###############################
### PREPROCESSING FUNCTIONS ###
###############################

def create_csv(path, file_name, results_type):
    file = os.path.join(path, file_name)
    f = open(file, 'w')
    
    if results_type == 'memory_log':
        f.write('mem_usage\n')
    elif results_type == 'file_metadata':
        f.write('tile_count,file_name,file_size,pixel_count,nodata_percentage\n')
    elif results_type == 'sequential_test':
        f.write('parameter,ex_time,peak_mem_usage\n')
    elif results_type == 'crop_test':
        f.write('tile_count,ex_time,peak_mem_usage\n')
    elif results_type == 'compute_test':
        f.write('tile_count,processes,parameter,ex_time,peak_mem_usage\n')
    elif results_type == 'single_compute_test':
        f.write('tile_count,file_name,parameter,ex_time,peak_mem_usage\n')
    elif results_type == 'single_compute_averages_test':
        f.write('tile_count,parameter,avg_ex_time,std_ex_time,avg_peak_mem_usage,std_peak_mem_usage\n')
    elif results_type == 'mosaic_test':
        f.write('tile_count,parameter,ex_time,peak_mem_usage\n')
    else:
        print('Results type provided invalid.')
    
    f.close()

    return file

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def merge_dataframes(csv_file1, csv_file2, filter_dictionary1, filter_dictionary2, merge_list, export_csv):
    # Read in and filter first CSV file
    data1 = pd.read_csv(csv_file1)
    df1 = pd.DataFrame(data1)

    for filter_key in filter_dictionary1:
        df1 = df1[df1[filter_key] == filter_dictionary1[filter_key]]

    # Read in and filter second CSV file
    data2 = pd.read_csv(csv_file2)
    df2 = pd.DataFrame(data2)

    for filter_key in filter_dictionary2:
        df2 = df2[df2[filter_key] == filter_dictionary2[filter_key]]

    # Merge on specified columns
    merged_df = pd.merge(df1, df2, on=merge_list, how="inner")

    # Export dataframe to CSV file
    merged_df.to_csv(export_csv, index=False)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def generate_elevation_data(roi, dataset):
    # Download elevation files for a specific region of interest within a given dataset
    gts.fetch_dem(shapefile=roi, dataset=dataset, save_to_txt=False, download=True)

    # Build mosaic from DEMs
    gts.build_mosaic(input_folder='dem_tiles', output_file='mosaic.tif', description='elevation')
    
    # Reproject the mosaic to Projected Coordinate System (PCS) EPSG:26918 - NAD83 UTM Zone 18N
    gts.reproject(input_file='mosaic.tif', output_file='elevation.tif', projection='EPSG:26918', cleanup=False)

#########################
### TESTING FUNCTIONS ###
#########################

def run_sequential_test(elevation_file, parameter_list):
    elevation_file = gts.determine_if_path(elevation_file) # Update path to elevation file if needed
    results_file = create_csv(os.getcwd(), 'sequential_test.csv', 'sequential_test') # Create results file

    for parameter in parameter_list:
        # Start trackers
        mem_log_file = create_csv(os.getcwd(), parameter+'_sequential_mem_log.csv', 'memory_log')
        start_memory_tracking('sequential-memory-tracking', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file)
        start_time = time.time()
    
        # Convert elevation to SDAT
        saga_elevation_file = os.path.join(os.path.dirname(elevation_file), "elevation.sdat")
        gts.convert_file_format(elevation_file, saga_elevation_file, "SAGA")
    
        # Compute parameter
        cmd_base = ["saga_cmd", "-c=1"]
        saga_parameter_file = os.path.join(os.path.dirname(elevation_file), parameter+".sgrd")
        
        if parameter == "slope": 
            cmd_full = cmd_base + ["ta_morphometry", "0", "-ELEVATION", saga_elevation_file.replace(".sdat",".sgrd"), "-SLOPE", saga_parameter_file]
        elif parameter == "aspect":
            cmd_full = cmd_base + ["ta_morphometry", "0", "-ELEVATION", saga_elevation_file.replace(".sdat",".sgrd"), "-ASPECT", saga_parameter_file]
        elif parameter == "profile_curvature":
            cmd_full = cmd_base + ["ta_morphometry", "0", "-ELEVATION", saga_elevation_file.replace(".sdat",".sgrd"), "-C_PROF", saga_parameter_file]
        elif parameter == "plan_curvature":
            cmd_full = cmd_base + ["ta_morphometry", "0", "-ELEVATION", saga_elevation_file.replace(".sdat",".sgrd"), "-C_PLAN", saga_parameter_file]
        elif parameter == "hillshade":
            cmd_full = cmd_base + ["ta_lighting", "0", "-ELEVATION", saga_elevation_file.replace(".sdat",".sgrd"), "-SHADE", saga_parameter_file]
        elif parameter == "convergence_index":
            cmd_full = cmd_base + ["ta_morphometry", "1", "-ELEVATION", saga_elevation_file.replace(".sdat",".sgrd"), "-RESULT", saga_parameter_file]

        gts.bash(cmd_full) # Run
    
        # Convert parameter to GeoTIFF
        parameter_file = os.path.join(os.path.dirname(elevation_file), parameter+".tif")
        gts.convert_file_format(saga_parameter_file.replace(".sgrd",".sdat"), parameter_file, "GTiff")
    
        # End trackers
        end_time = time.time()
        end_memory_tracking('sequential-memory-tracking')

        # Write results of computing
        write_results_to_csv(results_file, [parameter, end_time-start_time, compute_peak_memory_usage(mem_log_file)])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------


def run_crop_test(elevation_file, tile_counts):
    elevation_file = gts.determine_if_path(elevation_file) # Update path to elevation file if needed
    results_file = create_csv(os.getcwd(), 'crop_test.csv', 'crop_test') # Create results file
    metadata_file = create_csv(os.getcwd(), 'elevation_metadata.csv', 'file_metadata') # Create metadata results file

    # Store metadata for elevation file before split
    file_size = os.path.getsize(elevation_file)
    nodata_percentage, pixel_count = get_nodata_percentage(elevation_file, -999999.0)
    write_results_to_csv(metadata_file, [1, os.path.basename(elevation_file), file_size, pixel_count, nodata_percentage])
    
    for tile_count in tile_counts:
        # Create path to store cropped tiles
        tile_path = os.path.join(os.getcwd(), str(tile_count) + '_tiles')
        Path(tile_path).mkdir(parents=True, exist_ok=True)
        
        # Start trackers
        mem_log_file = create_csv(tile_path, 'crop_mem_log.csv', 'memory_log')
        start_memory_tracking('crop-memory-tracking', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file)
        start_time = time.time()
        
        # Begin cropping
        elevation_tile_path = os.path.join(tile_path,'elevation_tiles')
        gts.crop_into_tiles(input_file=elevation_file, output_folder=elevation_tile_path, num_tiles=tile_count)

        # End trackers
        end_time = time.time()
        end_memory_tracking('crop-memory-tracking')

        # Write results of cropping
        write_results_to_csv(results_file, [tile_count, end_time-start_time, compute_peak_memory_usage(mem_log_file)])

        # Get some metadata related to the elevation tiles
        elevation_files = sorted(glob.glob(os.path.join(elevation_tile_path, "*.tif")))
        for elev_file in elevation_files:
            file_name = os.path.basename(elev_file)
            file_size = os.path.getsize(elev_file)
            nodata_percentage, pixel_count = get_nodata_percentage(elev_file, -999999.0)

            # Write metadata to file
            write_results_to_csv(metadata_file, [tile_count, file_name, file_size, pixel_count, nodata_percentage])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_compute_test(elevation_tile_folders, parameter_list, num_processes=[1]):
    results_file = create_csv(os.getcwd(), 'compute_test.csv', 'compute_test') # Create results files
    
    for tile_folder in elevation_tile_folders:
        tile_folder = gts.determine_if_path(tile_folder) # Update path to elevation tile folder if needed
        elevation_files = sorted(glob.glob(os.path.join(tile_folder, "*.tif")))
        tile_count = len(elevation_files)
        
        for num_process in num_processes:
            for parameter in parameter_list:
                # Start trackers
                mem_log_file = create_csv(os.path.dirname(tile_folder), str(num_process)+'_'+parameter+'_mem_log.csv', 'memory_log')
                start_memory_tracking('compute-memory-tracking', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file)
                start_time = time.time()
        
                # Compute with GEOtiled
                gts.compute_geotiled(tile_folder, [parameter], num_processes=num_process)
        
                # End trackers
                end_time = time.time()
                end_memory_tracking('compute-memory-tracking')
        
                # Write results of computing
                write_results_to_csv(results_file, [tile_count, num_process, parameter, end_time-start_time, compute_peak_memory_usage(mem_log_file)])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_single_tile_compute_test(elevation_tile_folders, parameter_list):
    results_file = create_csv(os.getcwd(), 'single_compute_test.csv', 'single_compute_test') # Create results files
    averages_file = create_csv(os.getcwd(), 'single_compute_averages_test.csv', 'single_compute_averages_test') # Create averages files
    
    for tile_folder in elevation_tile_folders:
        tile_folder = gts.determine_if_path(tile_folder) # Update path to elevation tile folder if needed
        elevation_files = sorted(glob.glob(os.path.join(tile_folder, "*.tif")))
        tile_count = len(elevation_files)
        
        for elev_file in elevation_files:
            # Isolate copy of file in its own folder for single computation
            temp_folder = os.path.join(os.path.dirname(tile_folder), 'temp')
            Path(temp_folder).mkdir(parents=True, exist_ok=True)
            file_name = os.path.basename(elev_file)
            shutil.copyfile(elev_file, os.path.join(temp_folder, file_name))
            
            for parameter in parameter_list:
                # Start trackers
                mem_log_file = create_csv(os.path.dirname(tile_folder), str(file_name)+'_'+parameter+'_mem_log.csv', 'memory_log')
                start_memory_tracking('single-compute-memory-tracking', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file)
                start_time = time.time()
        
                # Compute with GEOtiled
                gts.compute_geotiled(temp_folder, [parameter], num_processes=1)
        
                # End trackers
                end_time = time.time()
                end_memory_tracking('single-compute-memory-tracking')
        
                # Write results of computing
                write_results_to_csv(results_file, [tile_count, file_name, parameter, end_time-start_time, compute_peak_memory_usage(mem_log_file)])

            # Delete temp copy
            shutil.rmtree(temp_folder)
            shutil.rmtree(os.path.join(os.path.dirname(temp_folder), "saga_temp"))

        # Compute and save averages for the single compute test to a seperate file
        data = pd.read_csv(results_file)
        df = pd.DataFrame(data)
        df = df[df['tile_count'] == tile_count]

        for parameter in parameter_list:
            ex_times = df[df['parameter'] == parameter]['ex_time']
            mem_usages = df[df['parameter'] == parameter]['peak_mem_usage']

            ex_time_avg = ex_times.mean()
            ex_time_std = ex_times.std()
            mem_usage_avg = mem_usages.mean()
            mem_usage_std = mem_usages.std()

            write_results_to_csv(averages_file, [tile_count, parameter, ex_time_avg, ex_time_std, mem_usage_avg, mem_usage_std])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_mosaic_test(input_folders, parameter_list):
    results_file = create_csv(os.getcwd(), 'mosaic_test.csv', 'mosaic_test') # Create results files

    for input_folder in input_folders:
        input_folder = gts.determine_if_path(input_folder)

        for parameter in parameter_list:
            parameter_folder = os.path.join(input_folder, parameter+"_tiles")
            input_files = sorted(glob.glob(os.path.join(parameter_folder, "*.tif")))
            tile_count = len(input_files)

            # Start trackers
            mem_log_file = create_csv(input_folder, parameter+'_mosaic_mem_log.csv', 'memory_log')
            start_memory_tracking('mosaic-tracking', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file)
            start_time = time.time()

            # Run mosaicking
            gts.build_mosaic_buffer(parameter_folder, os.path.join(input_folder, parameter+".tif"))

            # End trackers
            end_time = time.time()
            end_memory_tracking('mosaic-tracking')

            # Write results of computing
            write_results_to_csv(results_file, [tile_count, parameter, end_time-start_time, compute_peak_memory_usage(mem_log_file)])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def run_full_test(storage_path, roi, dataset, tile_counts, parameter_list, process_counts):
    # Set working directory
    gts.set_working_directory(storage_path)

    # Generate the elevation data
    generate_elevation_data(roi, dataset)

    # Run the sequential test
    run_sequential_test('elevation.tif', parameter_list)
    
    # Run cropping test    
    run_crop_test('elevation.tif', tile_counts)

    # Setup some file paths as other inputs to later functions
    tile_paths = []
    elevation_tile_paths = []
    for tile_count in tile_counts:
        tile_paths.append(str(tile_count)+'_tiles')
        elevation_tile_paths.append(str(tile_count)+'_tiles/elevation_tiles')

    # Run computing test
    run_compute_test(elevation_tile_paths, parameter_list, process_counts)

    # Run individual tile computing test
    run_single_tile_compute_test(elevation_tile_paths, parameter_list)

    # Run mosaicking test
    run_mosaic_test(tile_paths, parameter_list)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def calculate_total_run_times(sequential_csv, crop_csv, compute_csv, mosaic_csv, output_csv):
    # Write header of CSV file
    f = open(output_csv, 'w')
    f.write('tile_count,process_count,parameter,ex_time\n')
    f.close()
    
    # Read the CSVs into pandas dataframe
    data_seq = pd.read_csv(sequential_csv)
    df_seq = pd.DataFrame(data_seq)
    
    data_crop = pd.read_csv(crop_csv)
    df_crop = pd.DataFrame(data_crop)
    
    data_compute = pd.read_csv(compute_csv)
    df_compute = pd.DataFrame(data_compute)
    
    data_mosaic = pd.read_csv(mosaic_csv)
    df_mosaic = pd.DataFrame(data_mosaic)
    
    # Get unique variable values not associated with execution time
    process_counts = data_compute['processes'].unique()
    tile_counts = data_compute['tile_count'].unique()
    parameters = data_compute['parameter'].unique()
    
    # Read in sequential results first to new file
    df_seq = df_seq.reset_index()
    f = open(output_csv, 'a')
    for index, row in df_seq.iterrows():
        f.write('1,1,{},{}\n'.format(row['parameter'],row['ex_time']))
    f.close()
    
    # Get sum of execution times of GEOtiled components and write to file
    for tile_count in tile_counts:
        df_crop_time = df_crop[df_crop['tile_count'] == tile_count]['ex_time']
        df_compute_temp1 = df_compute[df_compute['tile_count'] == tile_count]
        df_mosaic_temp = df_mosaic[df_mosaic['tile_count'] == tile_count]
        for param in parameters:
            df_compute_temp2 = df_compute_temp1[df_compute_temp1['parameter'] == param]
            df_mosaic_time = df_mosaic_temp[df_mosaic_temp['parameter'] == param]['ex_time']
            for pc in process_counts:
                df_compute_time = df_compute_temp2[df_compute_temp2['processes'] == pc]['ex_time']
    
                total_ex_time = df_crop_time.iloc[0] + df_mosaic_time.iloc[0] + df_compute_time.iloc[0]
                f = open(output_csv, 'a')
                f.write('{},{},{},{}\n'.format(tile_count,pc,param,total_ex_time))
                f.close()

###############################
### VISUALIZATION FUNCTIONS ###
###############################

def create_single_bar_graph(csv_file, xdata, ydata, filter_dictionary, xlabel, ylabel, plot_title, yscale=1, plot_color='c'):
    # Read in CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    # Filter specific rows of CSV file from specified filter dictionary
    for filter_key in filter_dictionary:
        df = df[df[filter_key] == filter_dictionary[filter_key]]
    
    # Load in and scale x and y data
    x = df[xdata]
    y = df[ydata] / yscale

    # Create bar plot
    barWidth = 0.75
    bar = np.arange(len(x)) 
    plt.bar(bar, y, color=plot_color, width=barWidth, edgecolor='grey', label=ylabel) 

    # Add other labels
    plt.xticks([r for r in range(len(x))], x)
    plt.xlabel(xlabel, fontweight='bold', fontsize=15)
    plt.ylabel(ylabel, fontweight='bold', fontsize=15)
    plt.title(plot_title, fontweight='bold', fontsize=15)

    # Show Plot
    plt.show()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def create_multi_bar_graph(csv_file, xdata, ydata, ydata_filter, filter_list, filter_dictionary, xlabel, ylabel, plot_title, yscale=1, ydata_std=None):
    # Read in CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    # Filter specific rows of CSV file from specified filter dictionary
    for filter_key in filter_dictionary:
        df = df[df[filter_key] == filter_dictionary[filter_key]]

    # Load in x data
    x_data = df[xdata]

    # Set width of bar
    barWidth = 1 / (len(filter_list) + 1)

    # Set position of bar on X axis 
    bar_length = len(df[df[ydata_filter] == filter_list[0]][ydata])
    br = np.arange(bar_length)

    # Iterate through each filtered piece of ydata
    color_index = 0
    color_list = ['r','g','b','y','o','m']
    for filt in filter_list:
        filtered_ydata = df[df[ydata_filter] == filt][ydata] / yscale

        # Add the bar to the plot
        plt.bar(br, filtered_ydata, color=color_list[color_index], width=barWidth, edgecolor='grey', label=filt) 
        if ydata_std is not None:
            filtered_ydata_std = df[df[ydata_filter] == filt][ydata_std] / yscale
            plt.errorbar(br, filtered_ydata, yerr=filtered_ydata_std, fmt='k', linestyle='')
        color_index = 0 if color_index == 5 else color_index + 1

        # Update position of bar on x axis
        br = [x + barWidth for x in br]

    # Add labels
    plt.xlabel(xlabel, fontweight='bold', fontsize=15)
    plt.xticks([r + barWidth for r in range(bar_length)], x_data.unique())
    plt.ylabel(ylabel, fontweight='bold', fontsize=15)
    plt.title(plot_title, fontweight='bold', fontsize=15)

    # Show plot
    plt.legend()
    plt.show() 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def create_box_plot(csv_file, xdata, ydata, filter_dictionary, xlabel, ylabel, plot_title, yscale=1):
    # Read in CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    # Filter specific rows of CSV file from specified filter dictionary
    for filter_key in filter_dictionary:
        df = df[df[filter_key] == filter_dictionary[filter_key]]
    
    # Load in and scale x and y data
    x = df[xdata].unique()
    y = []
    for x_val in x:
        y.append(df[df[xdata] == x_val][ydata] / yscale)

    # Create plot
    plt.boxplot(y, labels=x)

    # Add labels
    plt.xlabel(xlabel, fontweight='bold', fontsize=15)
    plt.ylabel(ylabel, fontweight='bold', fontsize=15)
    plt.title(plot_title, fontweight='bold', fontsize=15)

    # Show plot
    plt.show()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def create_correlation_plot(csv_file, xdata, ydata, filter_dictionary, xlabel, ylabel, plot_title, yscale=1):
    # Read in CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    # Filter specific rows of CSV file from specified filter dictionary
    for filter_key in filter_dictionary:
        df = df[df[filter_key] == filter_dictionary[filter_key]]

    # Load in and scale x and y data
    x = df[xdata]
    y = df[ydata] / yscale
    
    # Create plot
    plt.scatter(x, y)

    # Add line fit
    m, b = np.polyfit(x, y, deg=1)
    plt.axline(xy1=(0, b), slope=m)

    # Add labels
    plt.xlabel(xlabel, fontweight='bold', fontsize=15)
    plt.ylabel(ylabel, fontweight='bold', fontsize=15)
    plt.title(plot_title, fontweight='bold', fontsize=15)

    # Show plot
    plt.show()