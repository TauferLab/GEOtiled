# This script is designed to be used with the Concurrency_Optimization.ipynb notebook. 
# It contains compute functions used by GEOtiled slightly modified to collection timing metrics.
# Last Updated: 12/06/2024
# Author: Gabriel Laboy (@glaboy-vol)

###############
### IMPORTS ###
###############

from pathlib import Path
from osgeo import gdal

import matplotlib.pyplot as plt
import crop_functions as cp
import pandas as pd
import numpy as np

import multiprocessing
import geotiled
import shutil
import glob
import os

#################
### CONSTANTS ###
#################

ALL_PARAMETERS = ["slope", "aspect", "hillshade"]

############################
### FUNCTION DEFINITIONS ###
############################

def set_file_description(file, description):
    """
    Sets a band description for a raster file.

    Uses the GDAL library to set a description to the first band of a 
    raster file.

    Parameters
    ----------
    file : str
        Path to file to update description for.
    description : str
        Description to be added to the first band of the raster file.
    """
    
    # Set description of data to parameter name
    dataset = gdal.Open(file)
    band = dataset.GetRasterBand(1)
    band.SetDescription(description)
    dataset = None

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def gdal_compute_parameters_parallel(input_file, items):
    """
    Computes a list of parameters using the GDAL library.

    Computes a set of user-specified parameters using algorithms implemented by the GDAL library.
    Meant to be used with a multiprocessing pool.

    Parameters
    ----------
    input_file : str
        Path to input file to compute parameters from.
    items : List[]
        A list that should only contain a list of parameters to compute.
    """

    # Extract parameters to compute
    param_list = items[0]

    # Compute each parameter
    for param in param_list:
        # Configure path to output file
        output_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), param+'_tiles', os.path.basename(input_file))
        
        # Configure compute options
        dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
        if param == 'aspect':
            dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])

        # Compute
        gdal.DEMProcessing(output_file, input_file, processing=param, options=dem_options)

        # Set description of data to parameter name
        set_file_description(output_file, param)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def saga_compute_parameters_parallel(input_file, items):  
    """
    Computes a list of parameters using the SAGA library.

    Computes a set of user-specified parameters using algorithms implemented by the SAGA library.
    The function also takes care of converting the input file from GeoTIFF to SGRD, and vice versa
    for the parameter files produced. Meant to be used with a multiprocessing pool.

    Parameters
    ----------
    input_file : str
        Path to input file to compute parameters from.
    items : List[]
        A list that should only contain a list of parameters to compute.
    """

    # Extract parameters to compute
    param_list = items[0]

    # Convert input file to SGRD
    geotiled.convert_file_format(input_file, input_file.replace('.tif','.sdat'), 'SAGA')

    # Compute each parameter
    if 'hillshade' in param_list:
        # Configure command and run
        output_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'hillshade_tiles', os.path.basename(input_file))
        cmd = ["saga_cmd", "-c=1", "ta_lighting", "0", "-ELEVATION", input_file.replace('.tif','.sgrd'), "-SHADE", output_file.replace('.tif','.sgrd')]
        geotiled.bash(cmd)

        # Convert output file to GeoTIFF
        geotiled.convert_file_format(output_file.replace('.tif','.sdat'), output_file, 'GTiff')

        # Set description of data to parameter name
        set_file_description(output_file, 'hillshade')

    if any(param in param_list for param in ["slope","aspect"]):
        # Configure command and run
        cmd = ["saga_cmd", "-c=1", "ta_morphometry", "0", "-ELEVATION", input_file.replace('.tif','.sgrd')]
        output_slope_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'slope_tiles', os.path.basename(input_file))
        output_aspect_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'aspect_tiles', os.path.basename(input_file))
        if "slope" in param_list:
            cmd = cmd + ["-SLOPE", output_slope_file.replace('.tif','.sgrd')]
        if "aspect" in param_list:
            cmd = cmd + ["-ASPECT", output_aspect_file.replace('.tif','.sgrd')]
        geotiled.bash(cmd)

        # Convert output files to GeoTIFF
        if "slope" in param_list:
            geotiled.convert_file_format(output_slope_file.replace('.tif','.sdat'), output_slope_file, 'GTiff')
            set_file_description(output_slope_file, 'slope')
        if "aspect" in param_list:
            geotiled.convert_file_format(output_aspect_file.replace('.tif','.sdat'), output_aspect_file, 'GTiff')
            set_file_description(output_aspect_file, 'aspect')
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_compute(input_folder, parameter_list, method='SAGA', cleanup=False, verbose=False):
    """
    Configures the multiprocessing pool for GEOtiled to begin computing terrain parameters.
    
    This function utilizes the multiprocessing library to allow for computation of parameters on different elevation files at the same time.
    The 'all' keyword can be passed to compute all parameters.
    
    Parameters
    ----------
    input_folder : str
        Name/path of the folder containing DEM elevation files to compute terrain parameters from.
    parameter_list : List[str]
        List containing codes for terrain parameters to compute. The 'all' keyword will compute all terrain parameters.
    method : string, optional
        Determines if parameters should be computed with GDAL API instead of SAGA API (default is SAGA).
    cleanup : bool, optional
        Determine if elevation files should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """

    # Ensure input folder exists
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return
    
    # Check for 'all' keyword in parameter list or validate parameter entries
    if (method == 'GDAL') and "all" in parameter_list:
        parameter_list = ['slope', 'aspect', 'hillshade']
    elif "all" in parameter_list:
        parameter_list = ALL_PARAMETERS
    else:
        parameter_list = [parameter.lower() for parameter in parameter_list]
        for parameter in parameter_list:
            if parameter not in ALL_PARAMETERS: 
                print(parameter, "is an invalid parameter. Terminating execution.")
                return

    # Create folders for each parameter
    for parameter in parameter_list:
        parameter_folder = os.path.join(os.path.dirname(input_path), parameter+"_tiles")
        Path(parameter_folder).mkdir(parents=True, exist_ok=True)
    
    # Get all input files
    input_files = glob.glob(os.path.join(input_path, "*.tif"))
    
    items = []
    for input_file in input_files:
        items.append((input_file, [parameter_list]))
    
    # Setup multi-processing pool and compute
    num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    if method == 'GDAL':
        pool.starmap(gdal_compute_parameters_parallel, items)
    else:
        pool.starmap(saga_compute_parameters_parallel, items)

    # Remove files used to compute parameters
    if cleanup is True:
        if verbose is True: print("Cleaning files...")
        shutil.rmtree(input_path)

    # Successful completion message
    if verbose is True: print("GEOtiled computation done!")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def gdal_compute_parameters(input_file, parameter_list):
    """
    Computes a list of parameters using the GDAL library.

    Computes a set of user-specified parameters using algorithms implemented by the GDAL library.

    Parameters
    ----------
    input_file : str
        Path to input file to compute parameters from.
    parameter_list : List[str]
        List of parameters to compute.
    """

    # Compute each parameter
    for parameter in parameter_list:
        # Configure path to output file
        output_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), parameter+'_tiles', os.path.basename(input_file))
        
        # Configure compute options
        dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
        if parameter == 'aspect':
            dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
    
        # Compute
        gdal.DEMProcessing(output_file, input_file, processing=parameter, options=dem_options)
    
        # Set description of data to parameter name
        set_file_description(output_file, parameter)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def saga_compute_parameters(input_file, parameter_list):
    """
    Computes a list of parameters using the SAGA library.

    Computes a set of user-specified parameters using algorithms implemented by the SAGA library.
    The function also takes care of converting the input file from GeoTIFF to SGRD, and vice versa
    for the parameter files produced.

    Parameters
    ----------
    input_file : str
        Path to input file to compute parameters from.
    parameter_list : List[str]
        List of parameters to compute.
    """

    # Convert input file to SGRD
    geotiled.convert_file_format(input_file, input_file.replace('.tif','.sdat'), 'SAGA')

    # Compute each parameter
    if 'hillshade' in parameter_list:
        # Configure command and run
        output_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'hillshade_tiles', os.path.basename(input_file))
        cmd = ["saga_cmd", "-c=1", "ta_lighting", "0", "-ELEVATION", input_file.replace('.tif','.sgrd'), "-SHADE", output_file.replace('.tif','.sgrd')]
        geotiled.bash(cmd)

        # Convert output file to GeoTIFF
        geotiled.convert_file_format(output_file.replace('.tif','.sdat'), output_file, 'GTiff')

        # Set description of data to parameter name
        set_file_description(output_file, 'hillshade')

    if any(param in parameter_list for param in ["slope","aspect"]):
        # Configure command and run
        cmd = ["saga_cmd", "-c=1", "ta_morphometry", "0", "-ELEVATION", input_file.replace('.tif','.sgrd')]
        output_slope_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'slope_tiles', os.path.basename(input_file))
        output_aspect_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'aspect_tiles', os.path.basename(input_file))
        if "slope" in parameter_list:
            cmd = cmd + ["-SLOPE", output_slope_file.replace('.tif','.sgrd')]
        if "aspect" in parameter_list:
            cmd = cmd + ["-ASPECT", output_aspect_file.replace('.tif','.sgrd')]
        geotiled.bash(cmd)

        # Convert output files to GeoTIFF
        if "slope" in parameter_list:
            geotiled.convert_file_format(output_slope_file.replace('.tif','.sdat'), output_slope_file, 'GTiff')
            set_file_description(output_slope_file, 'slope')
        if "aspect" in parameter_list:
            geotiled.convert_file_format(output_aspect_file.replace('.tif','.sdat'), output_aspect_file, 'GTiff')
            set_file_description(output_aspect_file, 'aspect')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_and_compute_parallel(window, items):
    """
    Performs cropping and computing for a single tile.

    Crops a single tile from a larger raster file and computes a list of specified
    parameters from it using a specified method. The computed files have the added
    buffer region removed after computation. This function is meant to be used with
    a Python multiprocessing pool.

    Parameters
    ----------
    window : List[int]
        Coordinates specifying what part of the original raster data to extract.
    items : List[]
        List of important variables passed from multiprocessing pool. Items include
        the original raster file, the name of the cropped file, the buffer size, the
        method to compute the paramters with, and the list of parameters to compute,
        respectively.
    """

    # Sort out all items
    input_file = items[0]
    tile_file = items[1]
    buffer = items[2]
    method = items[3]
    params = items[4]

    # Crop tiles from input file and specified buffer and window
    window[0] = window[0] - buffer
    window[1] = window[1] - buffer
    window[2] = window[2] + (2*buffer)
    window[3] = window[3] + (2*buffer)
    cp.crop_pixels(input_file, tile_file, window)

    # Compute parameters
    if method == 'GDAL':
        gdal_compute_parameters(tile_file, params)
    else:
        saga_compute_parameters(tile_file, params)

    # Crop off buffer region of computed files
    for param in params:
        # Set paths
        buffered_param_file = os.path.join(os.getcwd(),param+'_tiles',os.path.basename(tile_file))
        unbuffered_param_file = os.path.join(os.getcwd(),'unbuffered_'+param+'_tiles',os.path.basename(tile_file))

        # Get buffered file data and create new file with buffer pixels cropped off
        ds = gdal.Open(buffered_param_file, 0)
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        param_window = [buffer, buffer, cols-(buffer*2), rows-(buffer*2)]
        cp.crop_pixels(buffered_param_file, unbuffered_param_file, param_window)
        
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_crop_and_compute(input_file, column_length, row_length, parameter_list, method='SAGA', buffer=10, cleanup=False, verbose=False):
    """
    Handles cropping of an elevation file, computation of parameters, and cropping of parameter buffer regions all in one multiprocessing pool.

    Computes windows for each tile from the original elevation file and passes window information along with other relevant information
    to a single multiprocessing pool which will handle cropping the tile, computing all requested parameters from it, and removing the 
    buffer region from the computed parameter files.

    Parameters
    ----------
    input_file : str
        Name/path of the elevation file to compute terrain parameters from.
    column_length : int 
        Number of columns (in pixels) for each tile.
    row_length : int 
        Number of rows (in pixels) for each tile.
    parameter_list : List[str]
        List containing codes for terrain parameters to compute. The 'all' keyword will compute all terrain parameters.
    method : string, optional
        Determines if parameters should be computed with GDAL API instead of SAGA API (default is SAGA).
    buffer : int
        Specifies the buffer size - overlapping pixels that is included in the borders between two tiles (default is 10).
    cleanup : bool, optional
        Determine if elevation files should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """

    # Get full path to input file if needed
    input_path = geotiled.determine_if_path(input_file)
    
    # Create folders to store data intermediate data in
    input_tiles = os.path.join(os.getcwd(),'elevation_tiles')
    Path(input_tiles).mkdir(parents=True, exist_ok=True)
    for parameter in parameter_list:
        Path(os.path.join(os.getcwd(),parameter+'_tiles')).mkdir(parents=True, exist_ok=True)
        Path(os.path.join(os.getcwd(),'unbuffered_'+parameter+'_tiles')).mkdir(parents=True, exist_ok=True)

    # Get the total number of rows and columns of the input file
    ds = gdal.Open(input_path, 0)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    
    # Get windows and file names of to-be-cropped files
    tile_info = []
    tile_count = 0
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
            
            # Add window and other relevant variables to items
            tile_file = os.path.join(input_tiles, "tile_{0:04d}.tif".format(tile_count))
            tile_info.append([window, tile_file])
            tile_count += 1 

    # Store all required variables for multiprocessing into list
    items = []
    for tile in tile_info:
        items.append((tile[0],[input_path,tile[1],buffer,method,parameter_list]))

    # Setup multi-processing pool and compute
    num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(crop_and_compute_parallel, items)

    # Cleanup cropped input tiles if specified
    if cleanup is True:
        if verbose is True: print("Cleaning intermediary files...")
        shutil.rmtree(input_tiles)

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

def plot_concurrency_results(csv_file):
    """
    Plots results of concurrency test.

    Plots results of separated vs unified concurrnecy in a single bar graph.
    Format of the csv file from the Concurrency_Optimization.ipynb notebook should be followed to prevent error.

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

    # Get mosaic time means for each method and tile size
    concurreny_means, concurreny_std = {}, {}
    for method in df['method'].unique().tolist():
        concurreny_execution_times, concurreny_std_times = [], []
        for ts in df['tile_size'].unique().tolist():
            concurreny_time = df[(df['method'] == method) & (df['tile_size'] == ts)]['execution_time'].mean()
            std = df[(df['method'] == method) & (df['tile_size'] == ts)]['execution_time'].std()
            concurreny_execution_times.append(round(concurreny_time,2))
            concurreny_std_times.append(std)
        concurreny_means[method] = concurreny_execution_times
        concurreny_std[method] = concurreny_std_times
    
    # Convert tile sizes from pixels to kilometers
    dimensions = []
    for ts in df['tile_size'].unique().tolist():
        dimensions.append(convert_pixels_to_km(ts, 30))
    
    # Plot figure
    x = np.arange(len(dimensions))
    width = (1/3)
    
    fig, ax = plt.subplots(layout='constrained')
    rect1 = ax.bar(x + width/2, concurreny_means['separated'], width, yerr=concurreny_std['separated'], label='Separated', color='firebrick')
    rect2 = ax.bar(x + 1.5*width, concurreny_means['unified'], width, yerr=concurreny_std['unified'], label='Unified', color='royalblue')
    ax.bar_label(rect1, fontsize=9, padding=3)
    ax.bar_label(rect2, fontsize=9, padding=3)
    ax.set_ylabel('Execution Time (s)')
    ax.set_xlabel('Tile Size (km)')
    ax.set_title('Time of Different Concurrency Methods based on Tile Size')
    ax.set_xticks(x + width, dimensions)
    ax.legend(loc='upper right', ncols=2)
    ax.set_ylim(0, int(df['execution_time'].max()+30))
    
    plt.show()