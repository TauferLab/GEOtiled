# This script is designed to be used with the Crop_Optimization.ipynb notebook. 
# It contains functions for cropping used by GEOtiled slightly modified to collection timing metrics.
# Last Updated: 12/06/2024
# Author: Gabriel Laboy (@glaboy-vol)

###############
### IMPORTS ###
###############

from pathlib import Path
from osgeo import gdal

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import multiprocessing
import geotiled
import os

############################
### FUNCTION DEFINITIONS ###
############################

def crop_pixels(input_file, output_file, window):
    """
    Crops a raster file to a specific region given a specified window.
    
    This function uses GDAL functions to crop data based off pixel coordinates rather than geospatial coordinates.

    Parameters
    ----------
    input_file : str
        Name of input file to be cropped.
    output_file : str
        Name of cropped output file.
    window : List[int], Tuple(int)
        List or tuple of format [left_x, top_y, width, height] where left_x and top_y are pixel coordinates of the upper-left corner 
        of the cropping window, and width and height specify the dimensions of the cropping window in pixels.
    """
    
    # Set full file paths for input and output
    input_path = geotiled.determine_if_path(input_file)
    output_path = geotiled.determine_if_path(output_file)
    
    # Window to crop by [left_x, top_y, width, height]
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

    # Perform translation
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_crop(input_file, output_folder, column_length, row_length, buffer=10, verbose=False):
    """
    Splits a GeoTIFF file into smaller, square-sized tiles.
    
    This function divides a GeoTIFF file into a specified sized with added buffer regions to maintain accuracy since
    computation of parameters requires using neighboring pixels.

    Parameters 
    ----------
    input_file : str
        Name/path of the GeoTIFF file in the data directory to crop.
    output_folder : str
        Name/path of the folder in the data directory to store the cropped tiles.
    column_length : int 
        Number of columns (in pixels) for each tile.
    row_length : int 
        Number of rows (in pixels) for each tile.
    buffer : int
        Specifies the buffer size - overlapping pixels that is included in the borders between two tiles (default is 10).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Update path to input file
    input_path = geotiled.determine_if_path(input_file)

    # Ensure file to crop exists
    if geotiled.validate_path_exists(input_path) == -1: return
    
    # Update path to out folder and create it if not done so
    output_path = geotiled.determine_if_path(output_folder)
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Get the total number of rows and columns of the input file
    ds = gdal.Open(input_path, 0)
    cols = ds.RasterXSize
    rows = ds.RasterYSize

    tile_count = 0 # Track number of tiles cropped

    # Begin cropping process
    for i in range(0, rows, row_length):
        # If the next iteration were to exceed the number of rows of the file, crop it off to the correct length
        nrows = row_length
        if i + row_length > rows:
            nrows = rows - i

        for j in range(0, cols, column_length):
            # If the next iteration were to exceed the number of columns of the file, crop it off to the correct length
            ncols = column_length
            if j + column_length > cols:
                ncols = cols - j

            # Create path to new tile file and set initial crop window
            tile_file = os.path.join(output_path, "tile_{0:04d}.tif".format(tile_count))
            window = [j, i, ncols, nrows]

            # Set coords of upper left corner with the buffer included
            window[0] = window[0] - buffer
            window[1] = window[1] - buffer

            # Set coords of bottom right corner with the buffer included
            window[2] = window[2] + buffer*2
            window[3] = window[3] + buffer*2
            
            # Crop the new tile
            crop_pixels(input_path, tile_file, window)
            if verbose is True: print(os.path.basename(tile_file), "cropped.")
            tile_count += 1 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_pixels_parallel(window, items):
    """
    Crops a raster file to a specific region given a specified window.
    
    This function uses GDAL functions to crop data based off pixel coordinates rather than geospatial coordinates.
    It is meant to be specially used with a python multiprocessing pool.

    Parameters
    ----------
    window : List[int], Tuple(int)
        List or tuple of format [left_x, top_y, width, height] where left_x and top_y are pixel coordinates of the upper-left corner 
        of the cropping window, and width and height specify the dimensions of the cropping window in pixels.
    items : List[str]
        List containing input path to file to crop and output path to file to save cropped file to, respectively.
    """
    
    # Set options and perform crop
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
    gdal.Translate(items[1], items[0], options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_crop(input_file, output_folder, column_length, row_length, buffer=10, verbose=False):
    """
    Splits a GeoTIFF file into smaller files, or tiles.
    
    This function divides a GeoTIFF file into smaller, size-specified files with added buffer regions to maintain accuracy since
    computation of parameters requires using neighboring pixels. This function also uses a python multiprocessing pool to 
    concurrently crop numerous files for improved execution time.

    Parameters 
    ----------
    input_file : str
        Name/path of the GeoTIFF file in the data directory to crop.
    output_folder : str
        Name/path of the folder in the data directory to store the cropped tiles.
    column_length : int 
        Number of columns (in pixels) for each tile.
    row_length : int 
        Number of rows (in pixels) for each tile.
    num_processes : int, optional
        Number of concurrent processes to use when cropping files (default is 2).
    buffer : int, optional
        Specifies the buffer size - overlapping pixels that is included in the borders between two tiles (default is 10).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Update path to input file
    input_path = geotiled.determine_if_path(input_file)

    # Ensure file to crop exists
    if geotiled.validate_path_exists(input_path) == -1: return
    
    # Update path to out folder and create it if not done so
    output_path = geotiled.determine_if_path(output_folder)
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Get the total number of rows and columns of the input file
    ds = gdal.Open(input_path, 0)
    cols = ds.RasterXSize
    rows = ds.RasterYSize

    tile_count = 0 # Track number of tiles cropped

    # Begin cropping process
    items = []
    for i in range(0, rows, row_length):
        # If the next iteration were to exceed the number of rows of the file, crop it off to the correct length
        nrows = row_length
        if i + row_length > rows:
            nrows = rows - i

        for j in range(0, cols, column_length):
            # If the next iteration were to exceed the number of columns of the file, crop it off to the correct length
            ncols = column_length
            if j + column_length > cols:
                ncols = cols - j

            # Create path to new tile file and set initial crop window
            window = [j, i, ncols, nrows]

            # Set coords of upper left corner with the buffer included
            window[0] = window[0] - buffer
            window[1] = window[1] - buffer

            # Set coords of bottom right corner with the buffer included
            window[2] = window[2] + buffer*2
            window[3] = window[3] + buffer*2
            
            # Add window and other relevant variables to items
            tile_file = os.path.join(output_path, "tile_{0:04d}.tif".format(tile_count))
            items.append((window, [input_path,tile_file]))
            tile_count += 1 

    # Concurrently crop tiles
    num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(crop_pixels_parallel, items)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def convert_pixels_to_km(tile_size, resolution):
    """
    Converts a string specifying tile size from pixels to kilometers.
    
    Converts tile size from pixel to kilometer dimension. Resolution of data should be specified in meters.

    Parameters 
    ----------
    tile_size : str
        Tile size (in pixels) specified in the format widthxheight.
    resolution : int
        Resolution of the data in meters.
    """
    
    dim = tile_size.split('x')
    length = round((int(dim[0]) * resolution) / 1000)
    width = round((int(dim[1]) * resolution) / 1000)
    return str(length)+'x'+str(width)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_crop_results(csv_file):
    """
    Plots results of cropping test.

    Plots results of sequential vs parallel cropping in a single bar graph.
    Format of the csv file from the Crop_Optimization.ipynb notebook should be followed to prevent error.

    Parameters
    ----------
    csv_file : str
        CSV file to read data from to plot.
    """
    
    # Update path to csv file
    file_path = geotiled.determine_if_path(csv_file)

    # Ensure csv file exists
    if geotiled.validate_path_exists(file_path) == -1: return

    # Read in CSV into pandas dataframe
    data = pd.read_csv(file_path)
    df = pd.DataFrame(data)
    
    # Get crop time means and std for each method and tile size
    crop_means, crop_std = {}, {}
    for method in df['method'].unique().tolist():
        crop_execution_times, crop_std_times = [], []
        for ts in df['tile_size'].unique().tolist():
            crop_time = df[(df['method'] == method) & (df['tile_size'] == ts)]['execution_time'].mean()
            std = df[(df['method'] == method) & (df['tile_size'] == ts)]['execution_time'].std()
            crop_execution_times.append(round(crop_time,2))
            crop_std_times.append(std)
        crop_means[method] = crop_execution_times
        crop_std[method] = crop_std_times
    
    # Convert tile sizes from pixels to kilometers
    dimensions = []
    for ts in df['tile_size'].unique().tolist():
        dimensions.append(convert_pixels_to_km(ts, 30))

    # Plot figure
    x = np.arange(len(dimensions))
    width = (1/3)
    
    fig, ax = plt.subplots(layout='constrained')
    
    rect1 = ax.bar(x + width/2, crop_means['sequential'], width, yerr=crop_std['sequential'], label='Sequential', color='lightskyblue')
    rect2 = ax.bar(x + 1.5*width, crop_means['parallel'], width, yerr=crop_std['parallel'], label='Parallel', color='plum')
    ax.bar_label(rect1, fontsize=9, padding=3)
    ax.bar_label(rect2, fontsize=9, padding=3)
    
    ax.set_ylabel('Execution Time (s)')
    ax.set_xlabel('Tile Size (km)')
    ax.set_title('Time to Crop based on Tile Size')
    ax.set_xticks(x + width, dimensions)
    ax.legend(loc='upper right', ncols=2)
    ax.set_ylim(0, int(df['execution_time'].max()+10))
    
    plt.show()