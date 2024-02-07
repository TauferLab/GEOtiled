"""
::

     ____    ____    _____   ____            __                            
    /\  _`\ /\  _`\ /\  __`\/\  _`\         /\ \__                         
    \ \ \L\_\ \ \L\_\ \ \/\ \ \ \L\_\  __  _\ \ ,_\  _ __    __      ____  
     \ \ \L_L\ \  _\L\ \ \ \ \ \  _\L /\ \/'\\ \ \/ /\`'__\/'__`\   /',__\ 
      \ \ \/, \ \ \L\ \ \ \_\ \ \ \L\ \/>  </ \ \ \_\ \ \//\ \L\.\_/\__, `\ 
       \ \____/\ \____/\ \_____\ \____//\_/\_\ \ \__\\ \_\\ \__/.\_\/\____/
        \/___/  \/___/  \/_____/\/___/ \//\/_/  \/__/ \/_/ \/__/\/_/\/___/ 
                                                                           
                                                                           

Library of extra GEOTiled functions that are not needed for the normal workflow. Compiled by Jay Ashworth
v0.0.1
GCLab 2023

Derived from original work by: Camila Roa (@CamilaR20), Eric Vaughan (@VaughanEric), Andrew Mueller (@Andym1098), Sam Baumann (@sam-baumann), David Huang (@dhuang0212), Ben Klein (@robobenklein)

See GEOTiled docs for more information.
"""

import os
from osgeo import gdal, ogr  # Install in a conda env: https://anaconda.org/conda-forge/gdal
import numpy as np
import pandas as pd
import math
import concurrent.futures

# In Ubuntu: sudo apt-get install grass grass-doc
# pip install grass-session
from grass_session import Session
import grass.script as gscript
import tempfile

# Increased the size of GDAL’s input-output buffer cache to reduce the number of look-up operations
gdal.SetConfigOption("GDAL_CACHEMAX", "512")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_params(input_prefix, parameters):
    """
    Compute various topographic parameters using GDAL and GRASS GIS.
    ----------------------------------------------------------------

    This function computes a range of topographic parameters such as slope, aspect, and hillshading for a given Digital Elevation Model (DEM) using GDAL and GRASS GIS libraries.

    Required Parameters
    -------------------
    input_prefix : str
        Prefix path for the input DEM (elevation.tif) and the resulting parameter files.
        For instance, if input_prefix is "/path/to/dem/", then the elevation file should be 
        "/path/to/dem/elevation.tif" and the resulting slope will be at "/path/to/dem/slope.tif", etc.
    parameters : list of str
        List of strings specifying which topographic parameters to compute. Possible values are:
        'slope', 'aspect', 'hillshading', 'twi', 'plan_curvature', 'profile_curvature', 
        'convergence_index', 'valley_depth', 'ls_factor'.

    Outputs
    -------
    None
        Files are written to the `input_prefix` directory based on the requested parameters.

    Notes
    -----
    - GDAL is used for slope, aspect, and hillshading computations.
    - GRASS GIS is used for other parameters including 'twi', 'plan_curvature', 'profile_curvature', and so on.
    - The function creates a temporary GRASS GIS session for processing.
    - Assumes the input DEM is named 'elevation.tif' prefixed by `input_prefix`.

    Error states
    ------------
    - If an unsupported parameter is provided in the 'parameters' list, it will be ignored.
    """

    # Slope
    if 'slope' in parameters:
        dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'], callback=gdal.TermProgress_nocb)
        gdal.DEMProcessing(input_prefix + 'slope.tif', input_prefix + 'elevation.tif', processing='slope', options=dem_options)
    # Aspect
    if 'aspect' in parameters:
        dem_options = gdal.DEMProcessingOptions(zeroForFlat=True, format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'], callback=gdal.TermProgress_nocb)
        gdal.DEMProcessing(input_prefix + 'aspect.tif', input_prefix + 'elevation.tif', processing='aspect', options=dem_options)
    # Hillshading
    if 'hillshading' in parameters:
        dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'], callback=gdal.TermProgress_nocb)
        gdal.DEMProcessing(input_prefix + 'hillshading.tif', input_prefix + 'elevation.tif', processing='hillshade', options=dem_options)

    # Other parameters with GRASS GIS
    if any(param in parameters for param in ['twi', 'plan_curvature', 'profile_curvature']):
        # define where to process the data in the temporary grass-session
        tmpdir = tempfile.TemporaryDirectory()

        s = Session()
        s.open(gisdb=tmpdir.name, location='PERMANENT', create_opts=input_prefix + 'elevation.tif')
        creation_options = 'BIGTIFF=YES,COMPRESS=LZW,TILED=YES' # For GeoTIFF files

        # Load raster into GRASS without loading it into memory (else use r.import or r.in.gdal)
        gscript.run_command('r.external', input=input_prefix + 'elevation.tif', output='elevation', overwrite=True)
        # Set output folder for computed parameters
        gscript.run_command('r.external.out', directory=os.path.dirname(input_prefix), format="GTiff", option=creation_options)

        if 'twi' in parameters:
            gscript.run_command('r.topidx', input='elevation', output='twi.tif', overwrite=True)

        if 'plan_curvature' in parameters:
            gscript.run_command('r.slope.aspect', elevation='elevation', tcurvature='plan_curvature.tif', overwrite=True)

        if 'profile_curvature' in parameters:
            gscript.run_command('r.slope.aspect', elevation='elevation', pcurvature='profile_curvature.tif', overwrite=True)

        if 'convergence_index' in parameters:
            gscript.run_command('r.convergence', input='elevation', output='convergence_index.tif', overwrite=True)

        if 'valley_depth' in parameters:
            gscript.run_command('r.valley.bottom', input='elevation', mrvbf='valley_depth.tif', overwrite=True)

        if 'ls_factor' in parameters:
            gscript.run_command('r.watershed', input='elevation', length_slope='ls_factor.tif', overwrite=True)


        tmpdir.cleanup()
        s.close()
        
        # Slope and aspect with GRASS GIS (uses underlying GDAL implementation)
        #vgscript.run_command('r.slope.aspect', elevation='elevation', aspect='aspect.tif', slope='slope.tif', overwrite=True)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_params_concurrently(input_prefix, parameters):
    """
    Compute various topographic parameters concurrently using multiple processes.
    ------------------------------------------------------------------------------

    This function optimizes the performance of the `compute_params` function by concurrently computing 
    various topographic parameters. It utilizes Python's concurrent futures for parallel processing.

    Required Parameters
    -------------------
    input_prefix : str
        Prefix path for the input DEM (elevation.tif) and the resulting parameter files.
        E.g., if `input_prefix` is "/path/to/dem/", the elevation file is expected at 
        "/path/to/dem/elevation.tif", and the resulting slope at "/path/to/dem/slope.tif", etc.
    parameters : list of str
        List of strings specifying which topographic parameters to compute. Possible values include:
        'slope', 'aspect', 'hillshading', 'twi', 'plan_curvature', 'profile_curvature', 
        'convergence_index', 'valley_depth', 'ls_factor'.

    Outputs
    -------
    None
        Files are written to the `input_prefix` directory based on the requested parameters.

    Notes
    -----
    - Utilizes a process pool executor with up to 20 workers for parallel computations.
    - Invokes the `compute_params` function for each parameter in the list concurrently.

    Error states
    ------------
    - Unsupported parameters are ignored in the `compute_params` function.
    - Potential for resource contention: possible if multiple processes attempt simultaneous disk writes or read shared input files.
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        for param in parameters:
            executor.submit(compute_params, input_prefix, param)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# def compute_params_saga(file_prefix):
#     # Compute 15 terrain parameters with saga inside tmux session
#     # tmux new-session -d -s "myTempSession" bash ./terrestrial_parameters.sh ./mosaic/TN_WGS84_30m_
#     command = ['tmux', 'new-session', '-d', '-s', 'terrainParamsSession', 'bash ./terrestrial_parameters.sh ' + file_prefix]
#     # command = ['bash', './terrestrial_parameters.sh', name_file]
#     bash(command)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_coord(input_file, output_file, upper_left, lower_right):
    """
    Crop a raster file to a specific region using upper-left and lower-right coordinates.
    --------------------------------------------------------------------------------------

    This function uses the GDAL library to crop a raster file based on specified coordinates.

    Required Parameters
    -------------------
    input_file : str
        Path to the input raster file intended for cropping.
    output_file : str
        Destination path where the cropped raster file will be saved.
    upper_left : tuple of float
        (x, y) coordinates specifying the upper-left corner of the cropping window.
        Must be in the same projection as the input raster.
    lower_right : tuple of float
        (x, y) coordinates specifying the lower-right corner of the cropping window.
        Must be in the same projection as the input raster.

    Outputs
    -------
    None
        Generates a cropped raster file saved at the designated `output_file` location.

    Notes
    -----
    - The `upper_left` and `lower_right` coordinates define the bounding box for cropping.
    - Employs GDAL's Translate method with specific creation options for cropping.
    - For shapefiles, ensure they are unzipped. Using zipped files can lead to GDAL errors.

    Error states
    ------------
    - GDAL might raise an error if provided coordinates fall outside the input raster's bounds.
    """

    # upper_left = (x, y), lower_right = (x, y)
    # Coordinates must be in the same projection as the raster
    window = upper_left + lower_right
    translate_options = gdal.TranslateOptions(projWin=window, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'],
                                              callback=gdal.TermProgress_nocb)
    gdal.Translate(output_file, input_file, options=translate_options)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_region(input_file, shp_file, output_file):
    """
    Crop a raster file based on a region defined by a shapefile.
    ------------------------------------------------------------------

    This function uses the GDAL library to crop a raster file according to the boundaries 
    specified in a shapefile.

    Required Parameters
    -------------------
    input_file : str
        Path to the input raster file intended for cropping.
    shp_file : str
        Path to the shapefile that outlines the cropping region.
    output_file : str
        Destination path where the cropped raster file will be saved.

    Outputs
    -------
    None
        Produces a cropped raster file at the designated `output_file` location using boundaries 
        from the `shp_file`.

    Notes
    -----
    - Utilizes GDAL's Warp method, setting the `cutlineDSName` option and enabling `cropToCutline` 
      for shapefile-based cropping.

    Error states
    ------------
    - GDAL may generate an error if the shapefile's boundaries exceed the input raster's limits.
    - GDAL can also report errors if the provided shapefile is invalid or devoid of geometries.
    """
    warp_options = gdal.WarpOptions(cutlineDSName=shp_file, cropToCutline=True, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'],
                                    callback=gdal.TermProgress_nocb)
    warp = gdal.Warp(output_file, input_file, options=warp_options)
    warp = None

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def extract_raster(csv_file, raster_file, band_names):
    """
    Extract raster values corresponding to the coordinates specified in the CSV file.
    --------------------------------------------------------------------------------------------

    This function reads the x and y coordinates from a CSV file and extracts the raster values 
    corresponding to those coordinates. The extracted values are added to the CSV file as new columns.

    Required Parameters
    -------------------
    csv_file : str
        String representing the path to the CSV file containing 'x' and 'y' columns with coordinates.
    raster_file : str
        String representing the path to the raster file to extract values from.
    band_names : list of str
        List of strings specifying the column names for each band's extracted values.

    Outputs
    -------
    None
        Modifies the provided CSV file to include new columns with extracted raster values based on band_names.

    Notes
    -----
    - The CSV file must contain columns named 'x' and 'y' specifying the coordinates.
    - The order of band_names should correspond to the order of bands in the raster_file.
    - Ensure that GDAL and pandas are properly installed to utilize this function.

    Error States
    ------------
    - If the CSV file does not have 'x' and 'y' columns, a KeyError will occur.
    - If the specified coordinates in the CSV file are outside the bounds of the raster, incorrect or no values may be extracted.
    """
    # Extract values from raster corresponding to
    df = pd.read_csv(csv_file)

    ds = gdal.Open(raster_file, 0)
    gt = ds.GetGeoTransform()

    n_bands = ds.RasterCount
    bands = np.zeros((df.shape[0], n_bands))

    for i in range(df.shape[0]):
        px = int((df['x'][i] - gt[0]) / gt[1])
        py = int((df['y'][i] - gt[3]) / gt[5])

        for j in range(n_bands):
            band = ds.GetRasterBand(j + 1)
            val = band.ReadAsArray(px, py, 1, 1)
            bands[i, j] = val[0]
    ds = None

    for j in range(n_bands):
        df[band_names[j]] = bands[:, j]

    df.to_csv(csv_file, index=None)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_extent(shp_file):
    """
    Get the bounding box (extent) of a shapefile.
    ----------------------------------------------

    This function extracts the extent or bounding box of a shapefile. The extent is returned as 
    the upper left and lower right coordinates.

    Required Parameters
    -------------------
    shp_file : str
        String representing the path to the shapefile.

    Outputs
    -------
    tuple of tuple
        Returns two tuples, the first representing the upper left (x, y) coordinate and the second 
        representing the lower right (x, y) coordinate.

    Notes
    -----
    - Ensure that OGR is properly installed to utilize this function.

    Error States
    ------------
    - If the provided file is not a valid shapefile or cannot be read, OGR may raise an error.
    """
    ds = ogr.Open(shp_file)
    layer = ds.GetLayer()
    ext = layer.GetExtent()
    upper_left = (ext[0], ext[3])
    lower_right = (ext[1], ext[2])

    return upper_left, lower_right   

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# def get_projection(input_file, output_file):
#     cmd = ['gdalsrsinfo', '-o', 'wkt', input_file, '>', output_file]
#     bash(cmd)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# def shp2csv(input_file, output_file):
#     cmd = ['ogr2ogr', '-f', 'CSV', output_file, input_file, '-lco', 'GEOMETRY=AS_XY']
#     bash(cmd)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def tif2csv(raster_file, band_names=['elevation'], output_file='params.csv'):
    """
    Convert raster values from a TIF file into CSV format.
    ------------------------------------------------------

    This function reads raster values from a specified raster TIF file and exports them into a CSV format.
    The resulting CSV file will contain columns representing the x and y coordinates, followed by columns
    for each band of data in the raster.

    Required Parameters
    -------------------
    raster_file : str
        Path to the input raster TIF file to be converted.

    Optional Parameters
    -------------------
    band_names : list
        Names for each band in the raster. The order should correspond to the bands in the raster file. 
        Default is ['elevation'].
    output_file : str
        Path where the resultant CSV file will be saved. Default is 'params.csv'.

    Outputs
    -------
    None
        The function will generate a CSV file saved at the specified `output_file` path, containing the raster values 
        and their corresponding x and y coordinates.

    Notes
    -----
    - The x and y coordinates in the output CSV correspond to the center of each pixel.
    - NaN values in the CSV indicate that there's no data or missing values for a particular pixel.

    Error States
    ------------
    - If the provided raster file is not present, invalid, or cannot be read, GDAL may raise an error.
    - If the number of provided `band_names` does not match the number of bands in the raster, the resulting CSV 
      might contain columns without headers or may be missing some data.
    """
    ds = gdal.Open(raster_file, 0)
    xmin, res, _, ymax, _, _ = ds.GetGeoTransform()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    xstart = xmin + res / 2
    ystart = ymax - res / 2

    x = np.arange(xstart, xstart + xsize * res, res)
    y = np.arange(ystart, ystart - ysize * res, -res)
    x = np.tile(x[:xsize], ysize)
    y = np.repeat(y[:ysize], xsize)

    n_bands = ds.RasterCount
    bands = np.zeros((x.shape[0], n_bands))
    for k in range(1, n_bands + 1):
        band = ds.GetRasterBand(k)
        data = band.ReadAsArray()
        data = np.ma.array(data, mask=np.equal(data, band.GetNoDataValue()))
        data = data.filled(np.nan)
        bands[:, k-1] = data.flatten()

    column_names = ['x', 'y'] + band_names
    stack = np.column_stack((x, y, bands))
    df = pd.DataFrame(stack, columns=column_names)
    df.dropna(inplace=True)
    df.to_csv(output_file, index=None)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






