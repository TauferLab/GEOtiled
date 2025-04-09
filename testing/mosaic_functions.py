# This script is designed to be used with the Mosaic_Optimization.ipynb notebook. 
# It contains functions for mosaic used by GEOtiled slightly modified to collection timing metrics.
# Last Updated: 12/23/2024
# Author: Gabriel Laboy (@glaboy-vol)

###############
### IMPORTS ###
###############

from pathlib import Path
from osgeo import gdal

import matplotlib.pyplot as plt
import compute_functions as cm
import crop_functions as cp
import pandas as pd
import numpy as np

import multiprocessing
import geotiled
import shutil
import glob
import os

############################
### FUNCTION DEFINITIONS ###
############################

def mosaic_average(input_folder, output_file, cleanup=False, verbose=False):
    """
    Builds mosaic from multiple GeoTIFF files that were cropped with a buffer region.

    This function is similar to the `build_mosaic` function but handles mosaicking together GeoTIFF files that were split to includes buffer regions
    by averaging points in the buffer regions together when merging.

    Parameters
    ----------
    input_folder : str
        Name/path of folder in data directory where files to mosaic together are located.
    output_file : str
        Name/path of mosaicked file produced.
    cleanup : bool, optional
        Determine if files used for mosaicking should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """
    
    # Check to ensure input folder exists
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return
    
    if verbose is True: print("Mosaicking process started...")
    
    # Get files from input folder to merge together
    input_files = glob.glob(os.path.join(input_path, "*.tif"))

    # Build full path for VRT and output file
    output_path = geotiled.determine_if_path(output_file)
    vrt_file = os.path.join(os.getcwd(), 'merged.vrt')

    if verbose is True: print("Building VRT...")
    vrt = gdal.BuildVRT(vrt_file, input_files)
    vrt = None  # closes file

    with open(vrt_file, "r") as f:
        contents = f.read()

    if verbose is True: print("Averaging buffer values...")
    if "<NoDataValue>" in contents:
        nodata_value = contents[contents.index("<NoDataValue>") + len(
            "<NoDataValue>"): contents.index("</NoDataValue>")]  # To add averaging function
    else:
        nodata_value = 0

    code = '''band="1" subClass="VRTDerivedRasterBand">
  <PixelFunctionType>average</PixelFunctionType>
  <PixelFunctionLanguage>Python</PixelFunctionLanguage>
  <PixelFunctionCode><![CDATA[
import numpy as np

def average(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize,raster_ysize, buf_radius, gt, **kwargs):
    data = np.ma.array(in_ar, mask=np.equal(in_ar, {}))
    np.mean(data, axis=0, out=out_ar, dtype="float32")
    mask = np.all(data.mask,axis = 0)
    out_ar[mask] = {}
]]>
  </PixelFunctionCode>'''.format(nodata_value, nodata_value)

    sub1, sub2 = contents.split('band="1">', 1)
    contents = sub1 + code + sub2

    with open(vrt_file, "w") as f:
        f.write(contents)

    # Do translation to mosaicked file with bash function
    cmd = ["gdal_translate", "-co", "COMPRESS=LZW", "-co", "TILED=YES", "-co", 
           "BIGTIFF=YES", "--config", "GDAL_VRT_ENABLE_PYTHON", "YES", vrt_file, output_path]
    cm.bash(cmd)

    # Remove intermediary files used to build mosaic
    if cleanup is True:
        if verbose is True: print("Cleaning files...")
        shutil.rmtree(input_path)
    os.remove(vrt_file)

    if verbose is True: print("Mosaic process completed.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_buffer_regions(input_folder, output_folder, buffer=10):
    """
    Crops buffer pixels from a set of specified tiles.
    
    Crops buffer pixels off of all GeoTIFF tiles in a specified folder and produces a new folder
    of the cropped tiles. The cropping is done in parallel.

    Parameters
    ----------
    input_folder : str
        Name/path of the folder in data directory where files to mosaic together are located.
    output_file : str
        Name/path of mosaic file produced.
    buffer : int
        Size of buffer region of tiles in pixels.
    """

    # Check to ensure input folder exists
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return

    # Create output folder
    output_path = geotiled.determine_if_path(output_folder)
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Get all input files
    input_files = glob.glob(os.path.join(input_path, "*.tif"))
    
    items = []
    for input_file in input_files:
        file_name = os.path.basename(input_file)
        output_file = os.path.join(output_path, file_name)

        # Create window and crop new file without buffer pixels
        ds = gdal.Open(input_file, 0)
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        window = [buffer, buffer, cols-(buffer*2), rows-(buffer*2)]
        items.append((window, [input_file, output_file]))

    # Concurrently crop buffers from tiles
    num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(cp.crop_pixels_parallel, items)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def mosaic_concat(input_folder, output_file, buffer=10, description=None, cleanup=False, verbose=False):
    """
    Builds a mosaic out of multiple GeoTIFF files.
    
    This function creates a mosaic from a list of GeoTIFF files utilizing the buildVRT() function from the GDAL library.

    Parameters
    ----------
    input_folder : str
        Name/path of the folder in data directory where files to mosaic together are located.
    output_file : str
        Name/path of mosaic file produced.
    buffer : int
        Size of buffer region of tiles in pixels.
    description : str, optional
        Description to add to output raster band of mosaic file (default is None).
    cleanup : bool, optional
        Determines if files from `input_folder` should be deleted after mosaic is complete (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Create paths to files needed for computation
    vrt_path = os.path.join(os.getcwd(), 'merged.vrt')
    mosaic_path = geotiled.determine_if_path(output_file)
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return

    # Remove buffer from files
    if verbose is True: print("Cropping buffer region from input files...")
    crop_buffer_regions(input_path, "unbuffered_files")
    unbuffered_path = geotiled.determine_if_path("unbuffered_files")
    
    # Check for valid folder and put input files into list
    if verbose is True: print("Getting input files...")
    input_files = glob.glob(os.path.join(unbuffered_path, "*.tif"))

    # Build VRT (mosaic files together)
    if verbose is True: print("Constructing VRT...")
    vrt = gdal.BuildVRT(vrt_path, input_files)
    translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"])
    gdal.Translate(mosaic_path, vrt, options=translate_options)
    vrt = None  # close file

    # Update band description with name of terrain parameter
    if description is not None:
        if verbose is True: print("Updating band description...")
        dataset = gdal.Open(mosaic_path)
        band = dataset.GetRasterBand(1)
        band.SetDescription(description)
        dataset = None  # close file

    # Delete intermediary tiles used to build mosaic
    if cleanup is True:
        if verbose is True: print("Cleaning intermediary files...")
        shutil.rmtree(input_path)
        shutil.rmtree(unbuffered_path)
    os.remove(vrt_path)

    if verbose is True: print("Mosaic process complete.")

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

def plot_mosaic_results(csv_file, parameter):
    """
    Plots results of mosaic test.

    Plots results of averaging vs crop and concat mosaic in a single bar graph.
    Format of the csv file from the Mosaic_Optimization.ipynb notebook should be followed to prevent error.

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
    mosaic_means, mosaic_std = {}, {}
    for method in df['method'].unique().tolist():
        mosaic_execution_times, mosaic_std_times = [], []
        for ts in df['tile_size'].unique().tolist():
            mosaic_time = df[(df['method'] == method) & (df['tile_size'] == ts) & (df['parameter'] == parameter)]['execution_time'].mean()
            std = df[(df['method'] == method) & (df['tile_size'] == ts) & (df['parameter'] == parameter)]['execution_time'].std()
            mosaic_execution_times.append(round(mosaic_time,1))
            mosaic_std_times.append(std)
        mosaic_means[method] = mosaic_execution_times
        mosaic_std[method] = mosaic_std_times
    
    # Convert tile sizes from pixels to kilometers
    dimensions = []
    for ts in df['tile_size'].unique().tolist():
        dimensions.append(convert_pixels_to_km(ts, 30))
    
    # Plot figure
    x = np.arange(len(dimensions))
    width = (1/3)
    
    fig, ax = plt.subplots(layout='constrained')
    rect1 = ax.bar(x + width/2, mosaic_means['averaging'], width, yerr=mosaic_std['averaging'], label='Averaging', color='darkorange')
    rect2 = ax.bar(x + 1.5*width, mosaic_means['concatenation'], width, yerr=mosaic_std['concatenation'], label='Crop+Concat', color='seagreen')
    ax.bar_label(rect1, fontsize=9, padding=3)
    ax.bar_label(rect2, fontsize=9, padding=3)
    ax.set_ylabel('Execution Time (s)')
    ax.set_xlabel('Tile Size (km)')
    ax.set_title('Time to Mosaic based on Tile Size')
    ax.set_xticks(x + width, dimensions)
    ax.legend(loc='upper right', ncols=2)
    ax.set_ylim(0, int(df['execution_time'].max()+30))
    
    plt.show()