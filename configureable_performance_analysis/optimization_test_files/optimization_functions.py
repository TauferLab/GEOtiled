from pathlib import Path
from osgeo import gdal
import multiprocessing
import pandas as pd
import subprocess
import geotiled
import shutil
import glob
import os

gdal.UseExceptions()

#################
### FUNCTIONS ###
#################

def bash(argv):
    arg_seq = [str(arg) for arg in argv] # Convert all arguments in list into a string
    proc = subprocess.Popen(arg_seq, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() # Synchronize
    stdout, stderr = proc.communicate() # Get standard output and error of command

    # Print error message if execution returned error
    if proc.returncode != 0:
        raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
            ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_pixels(input_file, output_file, window):
    # Set full file paths for input and output
    input_path = geotiled.determine_if_path(input_file)
    output_path = geotiled.determine_if_path(output_file)
    
    # Window to crop by [left_x, top_y, width, height]
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

    # Perform translation
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_crop(input_file, output_folder, column_length, row_length, buffer=10):
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
            tile_count += 1 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_pixels_parallel(window, items):
    # Set options and perform crop
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
    gdal.Translate(items[1], items[0], options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_crop(input_file, output_folder, column_length, row_length, num_processes=None, buffer=10):
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
    if num_processes is None:
        num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(crop_pixels_parallel, items)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def sequential_compute(input_file, method, parameter):
    if method == 'GDAL':
        if parameter == 'aspect':
            dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            gdal.DEMProcessing(f"{parameter}.tif", input_file, processing=parameter, options=dem_options)
        else:
            dem_options = gdal.DEMProcessingOptions(format='GTiff', alg='ZevenbergenThorne', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            gdal.DEMProcessing(f"{parameter}.tif", input_file, processing=parameter, options=dem_options)
    elif method == 'SAGA':
        if parameter == 'slope':
            cmd = ['saga_cmd', '-c=1', 'ta_morphometry', '0', '-ELEVATION', input_file, '-SLOPE', 'slope.tif']
            bash(cmd)
        elif parameter == 'aspect':
            cmd = ['saga_cmd', '-c=1', 'ta_morphometry', '0', '-ELEVATION', input_file, '-ASPECT', 'aspect.tif']
            bash(cmd)
        elif parameter == 'hillshade':
            cmd = ['saga_cmd', '-c=1', 'ta_lighting', '0', '-ELEVATION', input_file, '-SHADE', 'hillshade.tif']
            bash(cmd)
    
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

    # Compute each parameter
    if 'hillshade' in param_list:
        # Configure command and run
        output_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'hillshade_tiles', os.path.basename(input_file))
        cmd = ["saga_cmd", "-c=1", "ta_lighting", "0", "-ELEVATION", input_file, "-SHADE", output_file]
        bash(cmd)

    if any(param in param_list for param in ["slope","aspect"]):
        # Configure command and run
        cmd = ["saga_cmd", "-c=1", "ta_morphometry", "0", "-ELEVATION", input_file]
        output_slope_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'slope_tiles', os.path.basename(input_file))
        output_aspect_file = os.path.join(os.path.dirname(os.path.dirname(input_file)), 'aspect_tiles', os.path.basename(input_file))
        if "slope" in param_list:
            cmd = cmd + ["-SLOPE", output_slope_file]
        if "aspect" in param_list:
            cmd = cmd + ["-ASPECT", output_aspect_file]
        bash(cmd)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def parallel_compute(input_folder, parameter_list, method='SAGA', num_processes=None, cleanup=False):
    # Ensure input folder exists
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return

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
    if num_processes is None:
        num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    if method == 'GDAL':
        pool.starmap(gdal_compute_parameters_parallel, items)
    else:
        pool.starmap(saga_compute_parameters_parallel, items)

    # Remove files used to compute parameters
    if cleanup is True:
        shutil.rmtree(input_path)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def mosaic_average(input_folder, output_file, cleanup=False):
    # Check to ensure input folder exists
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return
    
    # Get files from input folder to merge together
    input_files = glob.glob(os.path.join(input_path, "*.tif"))

    # Build full path for VRT and output file
    output_path = geotiled.determine_if_path(output_file)
    vrt_file = os.path.join(os.getcwd(), 'merged.vrt')

    vrt = gdal.BuildVRT(vrt_file, input_files)
    vrt = None  # closes file

    with open(vrt_file, "r") as f:
        contents = f.read()

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
    bash(cmd)

    # Remove intermediary files used to build mosaic
    if cleanup is True:
        shutil.rmtree(input_path)
    os.remove(vrt_file)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_buffer_regions(input_folder, output_folder, num_processes=None, buffer=10):
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
    if num_processes is None:
        num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(crop_pixels_parallel, items)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def mosaic_concat(input_folder, output_file, num_processes=None, buffer=10, description=None, cleanup=False):
    # Create paths to files needed for computation
    vrt_path = os.path.join(os.getcwd(), 'merged.vrt')
    mosaic_path = geotiled.determine_if_path(output_file)
    input_path = geotiled.determine_if_path(input_folder)
    if geotiled.validate_path_exists(input_path) == -1: return

    # Remove buffer from files
    crop_buffer_regions(input_path, "unbuffered_tiles", num_processes)
    unbuffered_path = geotiled.determine_if_path("unbuffered_tiles")
    
    # Check for valid folder and put input files into list
    input_files = glob.glob(os.path.join(unbuffered_path, "*.tif"))

    # Build VRT (mosaic files together)
    vrt = gdal.BuildVRT(vrt_path, input_files)
    translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"])
    gdal.Translate(mosaic_path, vrt, options=translate_options)
    vrt = None  # close file

    # Update band description with name of terrain parameter
    if description is not None:
        dataset = gdal.Open(mosaic_path)
        band = dataset.GetRasterBand(1)
        band.SetDescription(description)
        dataset = None  # close file

    # Delete intermediary tiles used to build mosaic
    if cleanup is True:
        shutil.rmtree(input_path)
        shutil.rmtree(unbuffered_path)
    os.remove(vrt_path)