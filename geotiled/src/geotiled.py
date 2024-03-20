'''
GEOtiled Refactored Library v1.0.0
GCLab 2024

Compiled by Jay Ashworth (@washwor1) and Gabriel Laboy (@glaboy-vol)

Derived from original work by: Camila Roa (@CamilaR20), Eric Vaughan (@VaughanEric), Andrew Mueller (@Andym1098), Sam Baumann (@sam-baumann), David Huang (@dhuang0212), and Ben Klein (@robobenklein)

Learn more about GEOtiled from the paper: https://dl.acm.org/doi/pdf/10.1145/3588195.3595941
'''

from osgeo import osr, ogr, gdal
from pathlib import Path
from tqdm import tqdm

import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import numpy as np

import concurrent.futures
import multiprocessing
import subprocess
import requests
import zipfile
import shutil
import math
import glob
import os

# To install in Ubuntu: 1) sudo apt-get install grass grass-doc 2) pip install grass-session
# Need to do some work to get extensions with scripts
# https://github.com/OSGeo/grass-addons (repo with addon code and installation)
# https://grass.osgeo.org/grass83/manuals/g.extension.html (how to use g.extension)
# https://grass.osgeo.org/grass83/manuals/addons/r.valley.bottom.html (about valley depth script)
# https://grasswiki.osgeo.org/wiki/GRASS-QGIS_relevant_module_list (default script list - look at scripts starting with 'r.')
# from grass_session import Session
# import grass.script as gscript
# import tempfile

# USGS dataset codes used for fetch_dem function
USGS_DATASET_CODES = {"30m":"National Elevation Dataset (NED) 1 arc-second Current",
                      "10m":"National Elevation Dataset (NED) 1/3 arc-second Current"}

# GRASS modules needed for the parameters collected from https://grasswiki.osgeo.org/wiki/Terrain_analysis
# Look at r.stream.* for channel network stuff https://grasswiki.osgeo.org/wiki/Hydrological_Sciences
# For relative slope position: https://grass.osgeo.org/grass83/manuals/addons/r.slope.direction.html
TERRAIN_PARAMETER_CODES = {"SLP":"slope", 
                           "ASP":"aspect", 
                           "HLSD":"hillshade"}#, 
                           # "CNBL":"channel_network_base_level",
                           # "CND":"channel_network_distance",
                           # "CD":"closed_depressions",
                           # "CI":"convergence_index",
                           # "LSF":"ls_factor",
                           # "PLC":"plan_curvature",
                           # "PFC":"profile_curvature", 
                           # "RSP":"relative_slope_position", 
                           # "TCA":"total_catchment_area", 
                           # "TWI":"topographic_wetness_index",
                           # "VD":"valley_depth"}

TERRAIN_PARAMETER_SCRIPT_CODES = {"slope":"slope", 
                                  "aspect":"aspect", 
                                  "hillshade":"hillshade", 
                                  "channel_network_base_level":"r.stream.channel", #unsure
                                  "channel_network_distance":"r.stream.distance", #unsure
                                  "closed_depressions":"r.fill.dir", # has two outputs https://grass.osgeo.org/grass83/manuals/r.fill.dir.html
                                  "convergence_index":"r.convergence",
                                  "ls_factor":"r.watershed",
                                  "plan_curvature":"r.slope.aspect",
                                  "profile_curvature":"r.slope.aspect", 
                                  "relative_slope_position":"r.slope.direction", #unsure
                                  "total_catchment_area":"r.catchment", #probably
                                  "topographic_wetness_index":"r.topidx",
                                  "valley_depth":"r.valley.bottom"}

REGION_CODES = {"AL":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Alabama_State_Shape.zip", 
                "AK":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Alaska_State_Shape.zip", 
                "AZ":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Arizona_State_Shape.zip", 
                "AR":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Arkansas_State_Shape.zip", 
                "CA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_California_State_Shape.zip", 
                "CO":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Colorado_State_Shape.zip", 
                "CT":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Connecticut_State_Shape.zip", 
                "DE":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Delaware_State_Shape.zip", 
                "DC":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_District_of_Columbia_State_Shape.zip",
                "FL":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Florida_State_Shape.zip", 
                "GA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Georgia_State_Shape.zip", 
                "GU":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Guam_State_Shape.zip",
                "HI":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Hawaii_State_Shape.zip", 
                "ID":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Idaho_State_Shape.zip",
                "IL":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Illinois_State_Shape.zip", 
                "IN":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Indiana_State_Shape.zip", 
                "IA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Iowa_State_Shape.zip", 
                "KS":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Kansas_State_Shape.zip", 
                "KY":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Kentucky_State_Shape.zip",
                "LA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Louisiana_State_Shape.zip", 
                "ME":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Maine_State_Shape.zip", 
                "MD":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Maryland_State_Shape.zip", 
                "MA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Massachusetts_State_Shape.zip", 
                "MI":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Michigan_State_Shape.zip",
                "MN":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Minnesota_State_Shape.zip", 
                "MS":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Mississippi_State_Shape.zip", 
                "MO":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Missouri_State_Shape.zip", 
                "MT":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Montana_State_Shape.zip", 
                "NE":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Nebraska_State_Shape.zip",
                "NV":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Nevada_State_Shape.zip", 
                "NH":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_Hampshire_State_Shape.zip", 
                "NJ":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_Jersey_State_Shape.zip", 
                "NM":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_Mexico_State_Shape.zip", 
                "NY":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_York_State_Shape.zip",
                "NC":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_North_Carolina_State_Shape.zip", 
                "ND":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_North_Dakota_State_Shape.zip", 
                "MP":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Commonwealth_of_the_Northern_Mariana_Islands_State_Shape.zip",
                "OH":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Ohio_State_Shape.zip", 
                "OK":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Oklahoma_State_Shape.zip",
                "OR":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Oregon_State_Shape.zip", 
                "PA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Pennsylvania_State_Shape.zip", 
                "PR":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Puerto_Rico_State_Shape.zip", 
                "RI":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Rhode_Island_State_Shape.zip", 
                "SC":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_South_Carolina_State_Shape.zip",
                "SD":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_South_Dakota_State_Shape.zip", 
                "TN":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Tennessee_State_Shape.zip", 
                "TX":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Texas_State_Shape.zip", 
                "UT":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Utah_State_Shape.zip", 
                "VT":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Vermont_State_Shape.zip",
                "VA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Virginia_State_Shape.zip", 
                "VI":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_United_States_Virgin_Islands_State_Shape.zip", 
                "WA":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Washington_State_Shape.zip", 
                "WV":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_West_Virginia_State_Shape.zip", 
                "WI":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Wisconsin_State_Shape.zip",
                "WY":"https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Wyoming_State_Shape.zip"}

# Path to store all data generated by functions
Data_Directory = './'

# Used to silence a deprecation warning. 
gdal.UseExceptions()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def bash(argv):
    '''
    Execute bash commands using Popen.
    ----------------------------------

    This function acts as a wrapper to execute bash commands using the subprocess Popen method. Commands are executed synchronously, 
    and errors are caught and raised.

    Required Parameters
    -------------------
    argv : List
        List of arguments for a bash command. They should be in the order that you would arrange them in the command line (e.g., ["ls", "-lh", "~/"]).

    Outputs
    -------
    None
        The function has no outputs, with the exception of any outputs produced by the bash command.

    Returns
    -------
    None
        The function does not return any value.

    Error States
    ------------
    RuntimeError
        Raises a RuntimeError if Popen returns with an error, detailing the error code, stdout, and stderr.

    Notes
    -----
    - It's essential to ensure that the arguments in the 'argv' list are correctly ordered and formatted for the desired bash command.
    '''
    arg_seq = [str(arg) for arg in argv]
    proc = subprocess.Popen(arg_seq, stdout=subprocess.PIPE, stderr=subprocess.PIPE)#, shell=True)
    proc.wait() #... unless intentionally asynchronous
    stdout, stderr = proc.communicate()

    # Error catching: https://stackoverflow.com/questions/5826427/can-a-python-script-execute-a-function-inside-a-bash-script
    if proc.returncode != 0:
        raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
            ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def set_data_directory(path):
    '''
    Set path where all data generated from functions will be stored.
    ----------------------------------

    This function sets a global that stores the root directory where all data generated by functions in this script will store data.

    Required Parameters
    -------------------
    path : str
        String containing path to store data in.

    Outputs
    -------
    None
        The function has no outputs.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - It's good practice to run this function before executing anything else in the workflow.
    - If not set, data will be stored in the working directory where functions are being executed.
    '''
    # Create the directory if it doesn't exist yet
    Path(path).mkdir(parents=True, exist_ok=True)
    
    # Set the global
    global Data_Directory
    Data_Directory = path

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
    
def get_file_size(url):
    '''
    Retrieve the size of a file at a given URL in bytes.
    ----------------------------------------------------

    This function sends a HEAD request to the provided URL and reads the 'Content-Length' header to determine the size of the file. 
    It's primarily designed to support the `download_files` function to calculate download sizes beforehand.

    Required Parameters
    -------------------
    url : str
        String representing the URL from which the file size needs to be determined.

    Outputs
    -------
    None
        The function has no outputs.

    Returns
    -------
    int
        Size of the file at the specified URL in bytes. Returns 0 if the size cannot be determined.

    Notes
    -----
    - This function relies on the server's response headers to determine the file size.
    - If the server doesn't provide a 'Content-Length' header or there's an error in the request, the function will return 0.
    - This function's primary use is with `download_files()`. 
    '''
    try:
        response = requests.head(url)
        return int(response.headers.get('Content-Length', 0))
    except:
        return 0

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_file(url, folder, pbar):
    '''
    Download a single file given its URL and store it in a specified path.
    -----------------------------------------------------------------------

    This is a utility function that facilitates the downloading of files, especially within iterative download operations.

    Required Parameters
    -------------------
    url : str
        String containing the URL of the file intended for downloading.
    folder : str
        String specifying the path where the downloaded file will be stored.
    pbar : tqdm object
        Reference to the tqdm progress bar, typically used in a parent function to indicate download progress.

    Outputs
    -------
    File
        Creates a file in the designated 'folder' upon successful download.

    Returns
    -------
    int
        Returns an integer representing the number of bytes downloaded.

    Notes
    -----
    - This function is meant to be used inside of `download_files()`. If you use it outside of that, YMMV. 
    - If the file already exists in the specified folder, no download occurs, and the function returns 0.
    - Utilizes the requests library for file retrieval and tqdm for progress visualization.
    '''
    local_filename = os.path.join(folder, url.split('/')[-1])
    if os.path.exists(local_filename):
        return 0

    response = requests.get(url, stream=True)
    downloaded_size = 0
    with open(local_filename, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
            downloaded_size += len(chunk)
            pbar.update(len(chunk))
            
    return downloaded_size

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_files(download_list, download_folder):
    '''
    Download one or multiple files from provided URLs and store them in a desired folder.
    ---------------------------------------------------

    This function allows for the simultaneous downloading of files using threading, and showcases download progress via a tqdm progress bar.

    Required Parameters
    -------------------
    download_list : str or list of str
        Can either be:
        1. A string specifying the name of a .txt file. This file should contain URLs separated by newlines.
        2. A list of strings where each string is a URL.
    download_folder : str
        String denoting the name of the folder where the downloaded files will be stored.

    Outputs
    -------
    Folder
        Creates and stores all downloaded files in the specified 'download_folder'.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - The function uses `ThreadPoolExecutor` from the `concurrent.futures` library to achieve multi-threaded downloads for efficiency.
    - The tqdm progress bar displays the download progress.
    - If the 'input' argument is a string, it's assumed to be the path to a .txt file containing URLs.
    - Will not download files if the file already exists, but the progress bar will not reflect it. 
    '''
    # Create directory to store downloaded files in
    download_folder = Data_Directory + download_folder
    Path(download_folder).mkdir(parents=True, exist_ok=True)

    # Update download list to be full path to file
    #download_list = Data_Directory + download_list + '.txt'

    # Determine if URLs are raw list or from a text file
    if type(download_list) is not list:
        download_list = Data_Directory + download_list + '.txt' # Update list to be full path to text file
        with open(download_list, 'r', encoding='utf8') as dsvfile:
            urls = [url.strip().replace("'$", "")
                    for url in dsvfile.readlines()]
    else:
        urls = download_list
    
    # Compute total size of files to download
    total_size = sum(get_file_size(url.strip()) for url in urls)

    # Create progress bar to track completion of downloads
    with tqdm(total=total_size, unit='B', unit_scale=True, ncols=100, desc='Downloading', colour='green') as pbar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            futures = [executor.submit(
                download_file, url, download_folder, pbar) for url in urls]
            for future in concurrent.futures.as_completed(futures):
                size = future.result()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_shape_files(codes):
    '''
    Download a shape file from a given code representing a region of choice.
    -------------------------------------------------------------------------

    This function downloads shape file from the TNM download site using their API. Files will be stored in a dedicated shape_file folder.

    Required Parameters
    -------------------
    codes : str
        String with comma-separated codes representing states or regions to download shape files for.

    Outputs
    -------
    Files
        The downloaded shape files stored and renamed for organization. 

    Returns
    -------
    None
        The function does not return a value.

    Notes
    -----
    - If the shape file arleady exists, it will simply return the path to it instead of redownloading it.

    Error States
    ------------
    - If the shape file is unable to get a successful response from TNM, a response error will be raised.
    '''
    # Ensure codes passed are valid codes
    for code in codes:
        if code not in REGION_CODES:
            print(code + ' is not a valid region code. Terminating execution.')
            return
    
    # Create path to shape file folder
    shape_path = Data_Directory + 'shape_files/'

    # Iterate through all passed code
    download_links = []
    download_codes = []
    for code in codes:
        # Create path to shape file
        shape_file_path = shape_path + code + '/' + code + '.shp'

        # Check if doesn't exist, add download URL to list
        if not os.path.isfile(shape_file_path):
            # Get download URL based off region code
            download_links.append(REGION_CODES[code])
            download_codes.append(code)

    # Do the downloads if the list contains URLs
    if len(download_links) > 0:
        download_files(download_links, 'shape_files')

        # Unzip and rename files
        iter = 0
        for link in download_links:
            # Isolate zip file name
            zip_name = os.path.basename(link)

            # Update path to zip file
            zip_path = shape_path + zip_name

            # Unzip and get correct region files
            original_file = 'GU_StateOrTerritory'
            with zipfile.ZipFile(zip_path, 'r') as file:
                file.extract('Shape/'+original_file+'.dbf', path=shape_path)
                file.extract('Shape/'+original_file+'.prj', path=shape_path)
                file.extract('Shape/'+original_file+'.shp', path=shape_path)
                file.extract('Shape/'+original_file+'.shx', path=shape_path)
            
            file.close()
            
            # Rename folder and files
            new_directory = shape_path + download_codes[iter] + '/'
            os.rename(shape_path + 'Shape/', new_directory)
            os.rename(new_directory+original_file+'.dbf', new_directory+download_codes[iter]+'.dbf')
            os.rename(new_directory+original_file+'.prj', new_directory+download_codes[iter]+'.prj')
            os.rename(new_directory+original_file+'.shp', new_directory+download_codes[iter]+'.shp')
            os.rename(new_directory+original_file+'.shx', new_directory+download_codes[iter]+'.shx')

            # Delete original zip file
            os.remove(zip_path)

            iter += 1

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_extent(shp_file):
    '''
    Get the bounding box (extent) of a shapefile.
    ----------------------------------------------

    This function extracts the extent or bounding box of a shapefile. The extent is returned as the upper left and lower right coordinates.

    Required Parameters
    -------------------
    shp_file : str
        String representing the path to the shapefile.

    Outputs
    -------
    None
        This function has no outputs.

    Returns
    -------
    tuple of tuple
        Returns two tuples, the first representing the upper left (x, y) coordinate and the second representing the lower right (x, y) coordinate.

    Notes
    -----
    - Ensure that OGR is properly installed to utilize this function.

    Error States
    ------------
    - If the provided file is not a valid shapefile or cannot be read, OGR may raise an error.
    '''
    ds = ogr.Open(shp_file)
    layer = ds.GetLayer()
    ext = layer.GetExtent()
    upper_left = (ext[0], ext[3])
    lower_right = (ext[1], ext[2])

    return upper_left, lower_right  

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def fetch_dem(shape_file=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset="30m", prod_format="GeoTIFF", txt_file='download_urls', save_to_txt=True, download=False):
    '''
    Queries the USGS API for DEM data given specified parameters and optionally extracts download URLs.
    ----------------------------------------------------------------------------------------------------

    The function targets the USGS National Map API, fetching Digital Elevation Model (DEM) data based on provided parameters. It can automatically download these files or save their URLs to a .txt file.

    Optional Parameters
    -------------------
    shape_file : str list
        List with single code of shapefile with which a bounding box will be generated. Overrides the 'bbox' parameter if set. Default is None.
    bbox : dict
        Dictionary containing bounding box coordinates to query. Consists of xmin, ymin, xmax, ymax. Default is {"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}.
    dataset : str
        Specifies the USGS dataset to target. Default is 30m which correlates to "National Elevation Dataset (NED) 1/3 arc-second Current".
    prod_format : str
        Desired file format for the downloads. Default is "GeoTIFF".
    txt_file : str
        Name to save the .txt file containing URLs. Default is "download_urls".
    save_to_txt : bool
        Flag to determine if URLs should be saved to a .txt file. Default is True.
    download : bool
        If set to True, triggers automatic file downloads. Default is False.

    Outputs
    -------
    File
        Depending on configurations, saves URLs to a .txt file.
    Folder
        Depending on configurations, creates folder called 'dem_tiles' with downloaded elevation files using the `download_files` function.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - If both 'bbox' and 'shape_file' are provided, 'bbox' will take precedence.
    - Uses the USGS National Map API for data fetching. Ensure the chosen dataset and product format are valid.
    '''    
    # Get coordinate extents if a shape file was specified
    if shape_file is not None:
        if len(shape_file) > 1:
            print('Only one region code allowed. Terminating execution.')
            return
        
        print('Reading in shape file...')
        
        # Download shape file if it doesn't already exist
        download_shape_files(shape_file)
        
        # Create path shape file path
        shape_file_path = Data_Directory + 'shape_files/' + shape_file[0] + '/' + shape_file[0] + '.shp'

        # Get extents of shape file
        coords = get_extent(shape_file_path)
        bbox['xmin'] = coords[0][0]
        bbox['ymax'] = coords[0][1]
        bbox['xmax'] = coords[1][0]
        bbox['ymin'] = coords[1][1]

    # Construct the query parameters
    print('Setting boundary extents...')
    params = {
        "bbox": f"{bbox['xmin']},{bbox['ymin']},{bbox['xmax']},{bbox['ymax']}",
        "datasets": USGS_DATASET_CODES[dataset],
        "prodFormats": prod_format
    }

    # Make a GET request
    print('Requesting data from USGS...')
    base_url = "https://tnmaccess.nationalmap.gov/api/v1/products"
    response = requests.get(base_url, params=params)

    # Check for a successful request
    if response.status_code != 200:
        raise Exception(
            f"Failed to fetch data. Status code: {response.status_code}")

    # Convert JSON response to Python dict
    data = response.json()

    # Extract download URLs
    download_urls = [item['downloadURL'] for item in data['items']]

    # Save URLs to text file 
    if save_to_txt is True:
        print('Saving URLs to text file...')
        txt_Path = Data_Directory + txt_file + '.txt' # Build full path to text file
        with open(txt_Path, "w") as file:
            for url in download_urls:
                file.write(f"{url}\n")

    # Download the files from the URLs
    if download is True:
        download_files(download_urls, "dem_tiles")

    print('Fetch process complete.')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_mosaic(input_folder, output_file, description, cleanup=False):
    '''
    Build a mosaic of geo-tiles using the GDAL library.
    ---------------------------------------------------

    This function creates a mosaic from a list of geo-tiles. 
    It is an integral part of the GEOTILED workflow and is frequently used for merging tiles into a single mosaic file.

    Required Parameters
    -------------------
    input_folder : str
        String specifying folder name where tiles to build mosaic from are located.
    output_file : str
        String specifying name of file for mosaic produced.
    description : str
        String specifying description to attach to the output raster band. 

    Optional Parameters
    -------------------
    cleanup : bool
        Boolean specifying if tiles used to create mosaic should be deleted after mosaic is complete. Default is False.

    Outputs
    -------
    File
        Generates a .tif file that is the created mosaic.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - Ensure the input geo-tiles are compatible for merging.
    - The function utilizes the GDAL library's capabilities to achieve the desired mosaic effect.
    '''
    # Create paths to files being created
    vrt_path = Data_Directory + 'merged.vrt'
    mosaic_path = Data_Directory + output_file + '.tif'
    
    # Check for valid folder and put input files into list
    print('Getting input files...')
    if not Path(Data_Directory + input_folder).exists():
        print('The folder ' + Data_Directory + input_folder + ' does not exist. Terminating execution.')
        return
    input_files = glob.glob(Data_Directory + input_folder + '/*.tif')

    # Build VRT (mosaic)
    print('Constructing VRT...')
    vrt = gdal.BuildVRT(vrt_path, input_files)
    translate_options = gdal.TranslateOptions(creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES', 'NUM_THREADS=ALL_CPUS'])
    gdal.Translate(mosaic_path, vrt, options=translate_options)
    vrt = None  # close file

    # Update band description with name of terrain parameter
    print('Updating band description...')
    dataset = gdal.Open(mosaic_path)
    band = dataset.GetRasterBand(1)
    band.SetDescription(description)
    dataset = None  # close file

    # Delete intermediary tiles used to build mosaic
    if cleanup is True:
        print('Cleaning intermediary files...')
        shutil.rmtree(Data_Directory + input_folder)
        os.remove(vrt_path)

    print('Mosaic process complete.')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def reproject(input_file, output_file, projection, cleanup=False):
    '''
    Reproject a geospatial raster dataset (GeoTIFF) using the GDAL library.
    -------------------------------------------------------------------------

    This function reprojects a given GeoTIFF raster dataset from its original coordinate system to a new specified projection. The result is saved as a new raster file. The projection can be provided in multiple formats, including standard EPSG codes or WKT format.

    Required Parameters
    -------------------
    input_file : str
        String representing the filename of the GeoTIFF to be reprojected.
    output_file : str
        String representing the filename of the reprojected raster.
    projection : str
        String indicating the desired target projection. Input can be a EPSG code (e.g., EPSG:4326) or the path to a .wkt file.

    Optional Parameters
    -------------------
    cleanup : bool
        Boolean to determine if original file that was reprojected should be deleted after completion.

    Outputs
    -------
    File
        Generates a reprojected GeoTIFF file with the name specified by 'output_file'.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - The function supports multi-threading for improved performance on multi-core machines.
    - The source raster data remains unchanged; only a new reprojected output file is generated.
    '''
    # Set full file paths for input and output
    if Data_Directory not in input_file:
        input_file = Data_Directory + input_file + '.tif'
        output_file = Data_Directory + output_file + '.tif'  

    # Ensure file to reproject exists
    if not os.path.isfile(input_file):
        print(os.path.basename(input_file) + ' does not exist. Terminating execution.' )
        return
    
    print('Reprojecting ' + os.path.basename(input_file) + '...')
    
    # Set warp options for reprojection and warp
    warp_options = gdal.WarpOptions(dstSRS=projection, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES', 'NUM_THREADS=ALL_CPUS'],
                                    multithread=True, warpOptions=['NUM_THREADS=ALL_CPUS'])
    warp = gdal.Warp(output_file, input_file, options=warp_options)
    warp = None  # Close file

    # Delete files used for reprojection
    if cleanup is True:
        print('Cleaning intermediary files...')
        os.remove(input_file)
        os.remove(input_file + '.aux.xml')

    print('Reprojection process complete.')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_pixels(input_file, output_file, window):
    '''
    Crop a raster file to a specific region using provided pixel window.
    ---------------------------------------------------------------------

    This function uses the GDAL library to perform the cropping operation based on pixel coordinates rather than geospatial coordinates.

    Required Parameters
    -------------------
    input_file : str
        String representing the path to the input raster file to be cropped.
    output_file : str
        String representing the path where the cropped raster file should be saved.
    window : list or tuple
        List or tuple specifying the window to crop by in the format [left_x, top_y, width, height].
        Here, left_x and top_y are the pixel coordinates of the upper-left corner of the cropping window,
        while width and height specify the dimensions of the cropping window in pixels.

    Outputs
    -------
    File
        A cropped raster file saved at the specified output_file path.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - The function uses GDAL's Translate method with the `srcWin` option to perform the pixel-based cropping.
    - Must ensure that GDAL is properly installed to utilize this function.

    Error States
    ------------
    - If the specified pixel window is outside the bounds of the input raster, an error might be raised by GDAL.
    '''
    # Window to crop by [left_x, top_y, width, height]
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
    gdal.Translate(output_file, input_file, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_into_tiles(input_file, output_folder, num_tiles, buffer=10):
    '''
    Splits a mosaic image into smaller, equally-sized tiles and saves them to a specified folder.
    ----------------------------------------------------------------------------------------------

    The function divides the mosaic into a specified number of tiles (num_tiles), taking care to adjust the dimensions of edge tiles and add a buffer to each tile. 

    Required Parameters 
    --------------------
    input_file : str
        The name of the GeoTIFF file to crop. 
    output_folder : str
        The name of the folder where the cropped tile images will be saved.
    n_tiles : int 
        The total number of tiles to produce. Should be a perfect square number.

    Optional Parameters 
    --------------------
    buffer : int
        Specifies the size of the buffer region between tiles in pixels. Default is 10.

    Outputs
    -------
    Folder
        Folder with cropped tiles are saved to the specified 'output_folder'.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    ------
        - The function will automatically create a buffer of overlapping pixels that is included in the borders between two tiles. This is customizable with the "buffer" kwarg. 
    '''
    # Update path to input file
    input_file = Data_Directory + input_file + '.tif'

    # Ensure file to crop exists
    if not os.path.isfile(input_file):
        print(os.path.basename(input_file) + ' does not exist. Terminating execution.' )
        return
    
    # Update path to out folder and create it if not done so
    output_folder = Data_Directory + output_folder
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Square root number of tiles to help get even number of rows and columns
    num_tiles = math.sqrt(num_tiles)

    # Split rows and columns of original file into even number of pixels for total number of tiles specified
    ds = gdal.Open(input_file, 0)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x_win_size = int(math.ceil(cols / num_tiles))
    y_win_size = int(math.ceil(rows / num_tiles))

    tile_count = 0 # Track number of tiles cropped

    for i in range(0, rows, y_win_size):
        if i + y_win_size < rows:
            nrows = y_win_size
        else:
            nrows = rows - i

        for j in range(0, cols, x_win_size):
            if j + x_win_size < cols:
                ncols = x_win_size
            else:
                ncols = cols - j

            tile_file = output_folder + '/tile_' + '{0:04d}'.format(tile_count) + '.tif'
            win = [j, i, ncols, nrows]

            # Upper left corner
            win[0] = max(0, win[0] - buffer)
            win[1] = max(0, win[1] - buffer)

            w = win[2] + 2 * buffer
            win[2] = w if win[0] + w < cols else cols - win[0]

            h = win[3] + 2 * buffer
            win[3] = h if win[1] + h < rows else rows - win[1]  

            crop_pixels(input_file, tile_file, win)
            print('tile_' + '{0:04d}'.format(tile_count) + '.tif cropped.')
            tile_count += 1

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_params(input_file, param_list):
    '''
    Generate terrain parameters for an elevation model.
    ----------------------------------------------------------------------------------------------

    This function uses the GDAL and GRASS libraries to compute terrain parameters like slope, aspect, and hillshading from a provided elevation model in .tif format.

    Required Parameters 
    --------------------
    input_file : str
        The path to the file to compute other parameters on. File should have elevation data.
    param_list : str
        Comma seperated list of strings with parameters to compute.

    Outputs
    -------
    Files
        Saves computed tiles to the appropiate folders.
        
    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - GDAL is used for slope, aspect, and hillshading computations.
    - GRASS GIS is used for other parameters including 'twi', 'plan_curvature', 'profile_curvature', and so on.
    - The function creates a temporary GRASS GIS session for processing.
    - The generated parameter files adopt the following GDAL creation options: 'COMPRESS=LZW', 'TILED=YES', and 'BIGTIFF=YES'.
    '''
    # Get name of specific file being computed on
    input_file_name = os.path.basename(input_file)

    # Compute all terrain parameters
    for param in param_list:
        # Set directory where computed tile will be stored
        path = Data_Directory + param + '_tiles/'
        Path(path).mkdir(parents=True, exist_ok=True)
        output_file = path + input_file_name
        
        # Set correct options and compute parameter
        if param in ['slope', 'aspect', 'hillshade']:
            if param == 'aspect':
                dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
                gdal.DEMProcessing(output_file, input_file, processing=param, options=dem_options)
            else:
                dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
                gdal.DEMProcessing(output_file, input_file, processing=param, options=dem_options)
        # else:
        #     # Define where to process the data in the temporary grass-session
        #     tmpdir = tempfile.TemporaryDirectory()

        #     # Create GRASS session
        #     s = Session()
        #     s.open(gisdb=tmpdir.name, location='PERMANENT', create_opts=input_file)
        #     creation_options = 'BIGTIFF=YES,COMPRESS=LZW,TILED=YES' # For GeoTIFF files

        #     # Load raster into GRASS without loading it into memory (else use r.import or r.in.gdal)
        #     gscript.run_command('r.external', input=input_file, output='elevation', overwrite=True)
            
        #     # Set output folder for computed parameters
        #     gscript.run_command('r.external.out', directory=path, format="GTiff", option=creation_options)

        #     # Compute parameter
        #     if param == 'slope':
        #         gscript.run_command('r.slope.aspect', elevation='elevation', slope=input_file_name, overwrite=True)
        #     elif param == 'aspect':
        #         gscript.run_command('r.slope.aspect', elevation='elevation', aspect=input_file_name, overwrite=True)
        #     elif param == 'topographic_wetness_index':
        #         gscript.run_command('r.topidx', input='elevation', output=input_file_name, overwrite=True)
        #     elif param == 'plan_curvature':
        #         gscript.run_command('r.slope.aspect', elevation='elevation', tcurvature=input_file_name, overwrite=True)
        #     elif param == 'profile_curvature':
        #         gscript.run_command('r.slope.aspect', elevation='elevation', pcurvature=input_file_name, overwrite=True)
        #     elif param == 'convergence_index':
        #         gscript.run_command('r.convergence', input='elevation', output=input_file_name, overwrite=True) #addon
        #     elif param == 'valley_depth':
        #         gscript.run_command('r.valley.bottom', input='elevation', mrvbf=input_file_name, overwrite=True) #addon
        #     elif param == 'ls_factor':
        #         gscript.run_command('r.watershed', input='elevation', length_slope=input_file_name, overwrite=True) # Threshold required

        #     # Slope and aspect with GRASS GIS (uses underlying GDAL implementation)
        #     #vgscript.run_command('r.slope.aspect', elevation='elevation', aspect='aspect.tif', slope='slope.tif', overwrite=True)
            
        #     # Cleanup
        #     tmpdir.cleanup()
        #     s.close()
    
        # Update band description and nodata value (for GRASS params)
        dataset = gdal.Open(output_file)
        band = dataset.GetRasterBand(1)
        band.SetDescription(param)
        #band.SetNoDataValue(-9999)
        dataset = None

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_geotiled(input_folder, param_list, num_procs, cleanup=False):
    '''
    Creates multiprocessing pool for computing terrain parameters.
    ---------------------------------------------------

    This function uses the Python multiprocessing library to create numerous python instances to compute terrain parameters for multiple input tiles simultaneously.
    
    Required Parameters
    --------------------
    input_folder : str
        String containing the name of the folder containing files to compute parameters with. The tiles should contain elevation data.
    param_list : str list
        List of string codes of parameters to compute. 'all' keyword can be used instead to denote computing all parameters.
    num_procs : int
        Integer specifying number of multiprocessing instances to compute with.

    Optional Parameters
    --------------------
    cleanup : bool
        Boolean specifying if folder used to compute other params should be deleted. Default is False.

    Outputs
    -------
    Files
        Generates terrain parameter files at the specified paths for slope, aspect, and hillshading.
    
    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - Note that for the num_procs variable, it is better practice to set this to a smaller number if the system computing does not have a lot of RAM.
    '''
    # Check to ensure input folder exists
    if not Path(Data_Directory + input_folder).exists():
        print('The folder ' + Data_Directory + input_folder + ' does not exist. Terminating execution.')
        return
    
    # Check if all params are to be computed
    if 'all' in param_list:
        # Ensure only the 'all' code entered
        if len(param_list) > 1:
            print("Can't use 'all' keyword with other parameters in list")
            return
        
        param_list.remove('all')
        for key in TERRAIN_PARAMETER_CODES:
            param_list.append(key)

    # Ensure all codes entered are valid codes and put full parameter name in different list
    params = []
    for param in param_list:
        if param not in TERRAIN_PARAMETER_CODES:
            print("Param code " + param + " is not a valid code")
            return
        params.append(TERRAIN_PARAMETER_CODES[param])
            
    
    # Get files from input folder to compute on 
    print('Getting input files...')
    input_files = sorted(glob.glob(Data_Directory + input_folder + '/*.tif'))

    # Configure a list with the input files and selected params to support computing with pool.starmap()
    # https://superfastpython.com/multiprocessing-pool-starmap/
    items = []
    for input_file in input_files:
        items.append((input_file, params))

    # Create multiprocessing pool based off number of tiles to compute and compute params
    #print(params, input_files)
    print('Starting computation of parameters...')
    pool = multiprocessing.Pool(processes=num_procs) 
    pool.starmap(compute_params, items)

    # Remove files used to compute params
    if cleanup is True:
        print('Cleaning files...')
        shutil.rmtree(Data_Directory + input_folder)

    print('GEOtiled computation completed.')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_mosaic_filtered(input_folder, output_file, cleanup=False):
    '''
    Build a mosaic of geo-tiles using the GDAL library with added logic to handle overlapping regions.
    ---------------------------------------------------------------------------------------------------

    This function creates a mosaic from a list of geo-tiles and introduces extra logic to handle averaging when regions overlap.
    The function is similar to 'build_mosaic' but provides additional capabilities to ensure the integrity of data in overlapping regions.

    Required Parameters
    -------------------
    input_files : list of str
        List of strings containing paths to the geo-tiles to be merged.
    output_file : str
        String representing the desired location and filename for the output mosaic.

    Optional Parameters
    -------------------
    cleanup : bool
        Boolean that specifies if intermediary tiles used to build mosaic should be deleted after the process is complete. Default is False.

    Outputs
    -------
    File
        Generates a .tif file representing the created mosaic with the name specified by 'output_file'.

    Returns
    -------
        The function does not return any value.

    Notes
    -----
    - The function makes use of the GDAL library's capabilities and introduces Python-based pixel functions to achieve the desired averaging effect.
    - The function is particularly useful when there are multiple sources of geo-data with possible overlapping regions,
      ensuring a smooth transition between tiles.
    - Overlapping regions in the mosaic are handled by averaging pixel values.
    '''
    # Check to ensure input folder exists
    if not Path(Data_Directory + input_folder).exists():
        print('The folder ' + Data_Directory + input_folder + ' does not exist. Terminating execution.')
        return
    
    print('Mosaicking started for ' + output_file + '.tif...')
    
    # Get files from input folder to merge together
    input_files = glob.glob(Data_Directory + input_folder + '/*.tif')

    # Build full path for VRT and output file
    output_file = Data_Directory + output_file + '.tif'
    vrt_file = Data_Directory + 'merged.vrt'

    print('Building VRT...')
    vrt = gdal.BuildVRT(vrt_file, input_files)
    vrt = None  # closes file

    with open(vrt_file, 'r') as f:
        contents = f.read()

    print('Averaging buffer values...')
    if '<NoDataValue>' in contents:
        nodata_value = contents[contents.index('<NoDataValue>') + len(
            '<NoDataValue>'): contents.index('</NoDataValue>')]  # To add averaging function
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

    with open(vrt_file, 'w') as f:
        f.write(contents)

    # Do translation to mosaicked file with bash function
    cmd = ['gdal_translate', '-co', 'COMPRESS=LZW', '-co', 'TILED=YES', '-co', 
           'BIGTIFF=YES', '--config', 'GDAL_VRT_ENABLE_PYTHON', 'YES', vrt_file, output_file]
    bash(cmd)

    # Remove intermediary files used to build mosaic
    if cleanup is True:
        print('Cleaning files...')
        shutil.rmtree(Data_Directory + input_folder)
        os.remove(vrt_file)

    print('Mosaic process completed.')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_to_valid_data(input_file, output_file, block_size=512):
    '''
    Crops a border region of NaN values from a GeoTIFF file. 
    -----------------------------------------------------------

    Using a blocking method, the function will scan through the GeoTIFF file to determine 
    the extent of valid data in order to crop of excess border NaN values.

    Required Parameters 
    --------------------
    input_file : str
        Path to the input GeoTIFF file.
    output_file : str
        Desired path for the output cropped GeoTIFF file.

    Optional Parameters
    --------------------
    block_size : int
        Specifies the size of blocks used in computing the extent. Default is 512. This means that a max 512x512 pixel area will be loaded into memory at any time. 

    Outputs
    -------
    File
        Saves a cropped GeoTIFF file to the specified output path.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    ------
    - block_size is used to minimize RAM usage with a blocking technique. Adjust to fit your performance needs.
    '''
    src_ds = gdal.Open(input_file, gdal.GA_ReadOnly)
    src_band = src_ds.GetRasterBand(1)
    
    no_data_value = src_band.GetNoDataValue()
    
    gt = src_ds.GetGeoTransform()
    
    # Initialize bounding box variables to None
    x_min, x_max, y_min, y_max = None, None, None, None

    for i in range(0, src_band.YSize, block_size):
        # Calculate block height to handle boundary conditions
        if i + block_size < src_band.YSize:
            actual_block_height = block_size
        else:
            actual_block_height = src_band.YSize - i

        for j in range(0, src_band.XSize, block_size):
            # Calculate block width to handle boundary conditions
            if j + block_size < src_band.XSize:
                actual_block_width = block_size
            else:
                actual_block_width = src_band.XSize - j

            block_data = src_band.ReadAsArray(j, i, actual_block_width, actual_block_height)
            
            rows, cols = np.where(block_data != no_data_value)
            
            if rows.size > 0 and cols.size > 0:
                if x_min is None or j + cols.min() < x_min:
                    x_min = j + cols.min()
                if x_max is None or j + cols.max() > x_max:
                    x_max = j + cols.max()
                if y_min is None or i + rows.min() < y_min:
                    y_min = i + rows.min()
                if y_max is None or i + rows.max() > y_max:
                    y_max = i + rows.max()

    # Convert pixel coordinates to georeferenced coordinates
    min_x = gt[0] + x_min * gt[1]
    max_x = gt[0] + (x_max + 1) * gt[1]
    min_y = gt[3] + (y_max + 1) * gt[5]
    max_y = gt[3] + y_min * gt[5]
    
    out_ds = gdal.Translate(output_file, src_ds, projWin=[min_x, max_y, max_x, min_y], projWinSRS='EPSG:4326')
    
    out_ds = None
    src_ds = None

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_region(input_file, shape_file, output_file):
    '''
    Crop a raster file based on a region defined by a shapefile.
    ------------------------------------------------------------------

    This function uses the GDAL library to crop a raster file according to the boundaries specified in a shapefile.

    Required Parameters
    -------------------
    input_file : str
        Path to the input raster file intended for cropping.
    shape_file : str
        Path to the shapefile that outlines the cropping region.
    output_file : str
        Destination path where the cropped raster file will be saved.

    Outputs
    -------
    File
        Produces a cropped raster file at the designated 'output_file' location using boundaries from the 'shp_file'.

    Returns
    -------
        The function does not return any value.

    Notes
    -----
    - Utilizes GDAL's Warp method, setting the 'cutlineDSName' option and enabling 'cropToCutline' for shapefile-based cropping.

    Error states
    ------------
    - GDAL may generate an error if the shapefile's boundaries exceed the input raster's limits.
    - GDAL can also report errors if the provided shapefile is invalid or devoid of geometries.
    '''
    
    warp_options = gdal.WarpOptions(cutlineDSName=shape_file, cropToCutline=True, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
    warp = gdal.Warp(output_file, input_file, options=warp_options)
    warp = None

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def generate_img(tif, cmap='inferno', dpi=150, downsample=1, verbose=False, clean=False, title=None, 
                 nancolor='green', ztype="Z", zunit=None, xyunit=None, reproject_gcs=False, 
                 shp_files=None, crop_shp=False, bordercolor="black", borderlinewidth=1.5, saveDir = None):
    '''
    Plot a GeoTIFF image using matplotlib.
    --------------------------------------

    This function is a powerful plotting tool for GEOTiff files that uses GDAL, OSR, numpy, matplotlib, and geopandas. 
    We have tried to create a simple interface for the end user where you can input a tif file and an informational image will be generated.
    If the default image is not suited for your needs or if any of the information is incorrect, there are a series of keyword arguments that allow for user customizability. 

    Several major features that are not enabled by default include:

    - Automatic PCS to GCS conversion using the ``reproject_gcs`` flag.
    - Automated cropping with a shapefile using the ``shp_file`` parameter in addition to the ``crop_shp``.
    - Downsampling in order to reduce computation time using the ``downsample`` flag. 
    - A verbose mode that will print additional spatial information about the geotiff file using the ``verbose`` flag.
    - A clean mode that will print an image of the geotiff with no other information using the ``clean`` flag. 

    Required Parameters
    --------------------
    tif : str
        Name of the GeoTIFF file to plot.

    Optional Parameters
    -------------------
    cmap : str
        Colormap used for visualization. Default is 'inferno'.
    dpi : int
        Resolution in dots per inch for the figure. Default is 150.
    downsample : int
        Factor to downsample the image by. Default is 10.
    verbose : bool
        If True, print geotransform and spatial reference details. Default is False.
    clean : bool
        If True, no extra data will be shown besides the plot image. Default is False.
    title : str
        Title for the plot. Default will display the projection name.
    nancolor : str
        Color to use for NaN values. Default is 'green'.
    ztype : str
        Data that is represented by the z-axis. Default is 'Z'.
    zunit : str
        Units for the data values (z-axis). Default is None and inferred from spatial reference.
    xyunit : str
        Units for the x and y axes. Default is None and inferred from spatial reference.
    reproject_gcs : bool
        Reproject a given raster from a projected coordinate system (PCS) into a geographic coordinate system (GCS).
    shp_files : str list
        Comma-seperated list of strings with shape file codes to use for cropping. Default is None.
    crop_shp : bool
        Flag to indicate if the shapefiles should be used for cropping. Default is False.
    bordercolor : str
        Color for the shapefile boundary. Default is "black".
    borderlinewidth : float
        Line width for the shapefile boundary. Default is 1.5.

    Outputs
    -------
    Image
        Displays image visualizing inputed GeoTIFF data with specified parameters.
    
    Returns
    -------
    raster_array: np.ndarray
        Returns the raster array that was used for visualization. 

    Notes
    -----
    - Alternative colormaps can be found in the `matplotlib documentation <https://matplotlib.org/stable/users/explain/colors/colormaps.html>`_.
    - Shapemap cannot be in a .zip form. GDAL will throw an error if you use a .zip. We recommend using .shp. It can also cause issues if you don't have the accompanying files with the .shp file. (.dbf, .prj, .shx).
    - Must be used with Jupyter Notebooks to display results properly. Will Implement a feature to save output to dir eventually. 
    - Using ``shp_file`` without setting ``crop_shp`` will allow you to plot the outline of the shapefile without actually cropping anything. 
    '''
    # Initial setup
    tif_dir_changed = False

    # Update full path to tif file
    tif = Data_Directory + tif + '.tif'

    # Ensure file to plot exists
    if not os.path.isfile(tif):
        print(os.path.basename(tif) + ' does not exist. Terminating execution.' )
        return
    
    # Reproject raster into geographic coordinate system if needed
    if reproject_gcs:
        base_dir = os.path.dirname(tif)
        new_path = os.path.join(base_dir, "vis.tif")
        reproject(tif, new_path, "EPSG:4326")
        if crop_shp is False:
            new_path_crop = os.path.join(base_dir, "vis_trim_crop.tif")
            print("Cropping NaN values...")
            crop_to_valid_data(new_path,new_path_crop)
            os.remove(new_path)
            tif = new_path_crop
        else:
            tif = new_path
        tif_dir_changed = True

    # Crop using shapefiles if needed
    shape_paths = []
    if crop_shp and shp_files:
        # Check if the list is not empty
        if not shp_files:
            print("Shapefile list is empty. Skipping shapefile cropping.")
        else:
            # Download shape files if any are missing
            download_shape_files(shp_files)

            # Update path to shape files
            for i in range(len(shp_files)):
                shape_paths.append(Data_Directory + 'shape_files/' + shp_files[i] + '/' + shp_files[i] + '.shp')

            #print(shape_paths)
            # Read each shapefile, clean any invalid geometries, and union them
            gdfs = [gpd.read_file(shp_file).buffer(0) for shp_file in shape_paths]
            combined_geom = gdfs[0].unary_union
            for gdf in gdfs[1:]:
                combined_geom = combined_geom.union(gdf.unary_union)
            
            combined_gdf = gpd.GeoDataFrame(geometry=[combined_geom], crs=gdfs[0].crs)
                
            # Save the combined shapefile temporarily for cropping
            temp_combined_shp = os.path.join(os.path.dirname(tif), "temp_combined.shp")
            combined_gdf.to_file(temp_combined_shp)
                
            print("Cropping with combined shapefiles...")
            base_dir = os.path.dirname(tif)
            new_path = os.path.join(base_dir, "crop.tif")
            crop_region(tif, temp_combined_shp, new_path)
            if tif_dir_changed:
                os.remove(tif)
            tif = new_path
            tif_dir_changed = True
            
            # Remove the temporary combined shapefile
            os.remove(temp_combined_shp)
            os.remove(os.path.join(base_dir, "temp_combined.cpg"))
            os.remove(os.path.join(base_dir, "temp_combined.dbf"))
            os.remove(os.path.join(base_dir, "temp_combined.prj"))
            os.remove(os.path.join(base_dir, "temp_combined.shx"))

    print("Reading in tif for visualization...")
    dataset = gdal.Open(tif)
    band = dataset.GetRasterBand(1)

    geotransform = dataset.GetGeoTransform()
    spatial_ref = osr.SpatialReference(wkt=dataset.GetProjection())

    # Extract spatial information about raster
    proj_name = spatial_ref.GetAttrValue('PROJECTION')
    proj_name = proj_name if proj_name else "GCS, No Projection"
    data_unit = zunit or spatial_ref.GetLinearUnitsName()
    coord_unit = xyunit or spatial_ref.GetAngularUnitsName()
    z_type = ztype if band.GetDescription() == '' else band.GetDescription()


    if verbose:
        print(f"Geotransform:\n{geotransform}\n\nSpatial Reference:\n{spatial_ref}\n\nDocumentation on spatial reference format: https://docs.ogc.org/is/18-010r11/18-010r11.pdf\n")

    raster_array = gdal.Warp('', tif, format='MEM', 
                             width=int(dataset.RasterXSize/downsample), 
                             height=int(dataset.RasterYSize/downsample)).ReadAsArray()

    # Mask nodata values
    raster_array = np.ma.array(raster_array, mask=np.equal(raster_array, band.GetNoDataValue()))

    print("Plotting data...")

    # Set up plotting environment
    cmap_instance = plt.get_cmap(cmap)
    cmap_instance.set_bad(color=nancolor)

    # Determine extent
    ulx, xres, _, uly, _, yres = geotransform
    lrx = ulx + (dataset.RasterXSize * xres)
    lry = uly + (dataset.RasterYSize * yres)

    # Plot
    fig, ax = plt.subplots(dpi=dpi)
    sm = ax.imshow(raster_array, cmap=cmap_instance, vmin=np.nanmin(raster_array), vmax=np.nanmax(raster_array), 
                   extent=[ulx, lrx, lry, uly])
    if clean:
        ax.axis('off')
    else:
        # Adjust colorbar and title
        cbar = fig.colorbar(sm, fraction=0.046*raster_array.shape[0]/raster_array.shape[1], pad=0.04)
        cbar_ticks = np.linspace(np.nanmin(raster_array), np.nanmax(raster_array), 8)
        cbar.set_ticks(cbar_ticks)
        cbar.set_label(f"{z_type} ({data_unit}s)")

        ax.set_title(title if title else f"Visualization of GEOTiff data using {proj_name}.", fontweight='bold')
        ax.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, color='black', length=5, width=1)

        ax.set_title(title or f"Visualization of GEOTiff data using {proj_name}.", fontweight='bold')

    # Set up the ticks for x and y axis
    x_ticks = np.linspace(ulx, lrx, 5)
    y_ticks = np.linspace(lry, uly, 5)

    # Format the tick labels to two decimal places
    x_tick_labels = [f'{tick:.2f}' for tick in x_ticks]
    y_tick_labels = [f'{tick:.2f}' for tick in y_ticks]

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_yticklabels(y_tick_labels)

    # Determine x and y labels based on whether data is lat-long or projected
    y_label = f"Latitude ({coord_unit}s)" if spatial_ref.EPSGTreatsAsLatLong() else f"Northing ({coord_unit}s)"
    x_label = f"Longitude ({coord_unit}s)" if spatial_ref.EPSGTreatsAsLatLong() else f"Easting ({coord_unit}s)"
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)

    # ax.ticklabel_format(style='plain', axis='both')  # Prevent scientific notation on tick labels
    ax.set_aspect('equal')

    if shp_files:
        for shp_file in shape_paths:
            overlay = gpd.read_file(shp_file)
            overlay.boundary.plot(color=bordercolor, linewidth=borderlinewidth, ax=ax)

    print("Done. Image should appear soon...")

    if saveDir is not None:
        fig.savefig(saveDir)

    if tif_dir_changed:
        os.remove(tif)

    return raster_array

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_coord(input_file, output_file, upper_left, lower_right):
    '''
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
    File
        Generates a cropped raster file saved at the designated 'output_file' location.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - The `upper_left` and `lower_right` coordinates define the bounding box for cropping.
    - Employs GDAL's Translate method with specific creation options for cropping.
    - For shapefiles, ensure they are unzipped. Using zipped files can lead to GDAL errors.

    Error states
    ------------
    - GDAL might raise an error if provided coordinates fall outside the input raster's bounds.
    '''
    # upper_left = (x, y), lower_right = (x, y)
    # Coordinates must be in the same projection as the raster
    window = upper_left + lower_right
    translate_options = gdal.TranslateOptions(projWin=window, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])#,callback=gdal.TermProgress_nocb)
    gdal.Translate(output_file, input_file, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def change_raster_format(input_file, output_file, raster_format):
    '''
    Convert the format of a given raster file.
    -------------------------------------------

    This function leverages the GDAL library to facilitate raster format conversion.
    It allows users to convert the format of their .tif raster files to several supported formats, 
    specifically highlighting the GTiff and NC4C formats.

    Required Parameters
    -------------------
    input_file : str
        String containing the path to the input .tif file.
    output_file : str
        String representing the desired location and filename for the output raster.
    raster_format : str
        String indicating the desired output format for the raster conversion. 
        Supported formats can be found at GDAL Raster Formats: https://gdal.org/drivers/raster/index.html
        This function explicitly supports GTiff and NC4C.

    Outputs
    -------
    File
        A raster file in the desired format is generated at the location specified by 'output_file'.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - While GTiff and NC4C formats have been explicitly mentioned, 
      the function supports several other formats as listed in the GDAL documentation.
    - The function sets specific creation options for certain formats. 
      For example, the GTiff format will use LZW compression, tiling, and support for large files.
    '''
    # Supported formats: https://gdal.org/drivers/raster/index.html
    # SAGA, GTiff
    if raster_format == 'GTiff':
        translate_options = gdal.TranslateOptions(format=raster_format, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])#,callback=gdal.TermProgress_nocb)
    elif raster_format == 'NC4C':
        translate_options = gdal.TranslateOptions(format=raster_format, creationOptions=['COMPRESS=DEFLATE'])#,callback=gdal.TermProgress_nocb)
    else:
        translate_options = gdal.TranslateOptions(format=raster_format)#,callback=gdal.TermProgress_nocb)
    
    gdal.Translate(output_file, input_file, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def extract_raster(csv_file, raster_file, band_names):
    '''
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
    CSV
        Modifies the provided CSV file to include new columns with extracted raster values based on band_names.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - The CSV file must contain columns named 'x' and 'y' specifying the coordinates.
    - The order of band_names should correspond to the order of bands in the raster_file.
    - Ensure that GDAL and pandas are properly installed to utilize this function.

    Error States
    ------------
    - If the CSV file does not have 'x' and 'y' columns, a KeyError will occur.
    - If the specified coordinates in the CSV file are outside the bounds of the raster, incorrect or no values may be extracted.
    '''
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

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_stack(input_files, output_file):
    '''
    Stack a list of .tif files into a single .tif file with multiple bands.
    ------------------------------------------------------------------------

    This function takes multiple .tif files and combines them into a single .tif file where each input file represents a separate band. 
    This operation is useful when multiple datasets, each represented by a separate .tif file, need to be combined into a single multi-band raster.

    Required Parameters
    -------------------
    input_files : list of str
        List of strings containing paths to the .tif files to be stacked.
    output_file : str
        String representing the desired location and filename for the output stacked raster.

    Outputs
    -------
    File
        A multi-band .tif file is generated at the location specified by 'output_file'.

    Returns
    -------
        The function does not return any value.

    Notes
    -----
    - The function makes use of the GDAL library's capabilities to achieve the stacking operation.
    - Each input .tif file becomes a separate band in the output .tif file, retaining the order of the `input_files` list.
    '''
    # input_files: list of .tif files to stack
    vrt_options = gdal.BuildVRTOptions(separate=True)
    vrt = gdal.BuildVRT("stack.vrt", input_files, options=vrt_options)
    translate_options = gdal.TranslateOptions(creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])#,callback=gdal.TermProgress_nocb)
    gdal.Translate(output_file, vrt, options=translate_options)
    vrt = None  # closes file
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def tif2csv(raster_file, band_names=['elevation'], output_file='params.csv'):
    '''
    Convert raster values from a TIF file into CSV format.
    ------------------------------------------------------

    This function reads raster values from a specified raster TIF file and exports them into a CSV format.
    The resulting CSV file will contain columns representing the x and y coordinates, followed by columns for each band of data in the raster.

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
    CSV
        The function will generate a CSV file saved at the specified `output_file` path, containing the raster values and their corresponding x and y coordinates.

    Returns
    -------
        The function does not return any value.

    Notes
    -----
    - The x and y coordinates in the output CSV correspond to the center of each pixel.
    - NaN values in the CSV indicate that there's no data or missing values for a particular pixel.

    Error States
    ------------
    - If the provided raster file is not present, invalid, or cannot be read, GDAL may raise an error.
    - If the number of provided `band_names` does not match the number of bands in the raster, the resulting CSV 
      might contain columns without headers or may be missing some data.
    '''
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

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# def compute_params_concurrently(input_prefix, parameters):
#     """
#     Compute various topographic parameters concurrently using multiple processes.
#     ------------------------------------------------------------------------------

#     This function optimizes the performance of the `compute_params` function by concurrently computing 
#     various topographic parameters. It utilizes Python's concurrent futures for parallel processing.

#     Required Parameters
#     -------------------
#     input_prefix : str
#         Prefix path for the input DEM (elevation.tif) and the resulting parameter files.
#         E.g., if `input_prefix` is "/path/to/dem/", the elevation file is expected at 
#         "/path/to/dem/elevation.tif", and the resulting slope at "/path/to/dem/slope.tif", etc.
#     parameters : list of str
#         List of strings specifying which topographic parameters to compute. Possible values include:
#         'slope', 'aspect', 'hillshading', 'twi', 'plan_curvature', 'profile_curvature', 
#         'convergence_index', 'valley_depth', 'ls_factor'.

#     Outputs
#     -------
#     None
#         Files are written to the `input_prefix` directory based on the requested parameters.

#     Notes
#     -----
#     - Utilizes a process pool executor with up to 20 workers for parallel computations.
#     - Invokes the `compute_params` function for each parameter in the list concurrently.

#     Error states
#     ------------
#     - Unsupported parameters are ignored in the `compute_params` function.
#     - Potential for resource contention: possible if multiple processes attempt simultaneous disk writes or read shared input files.
#     """
#     with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
#         for param in parameters:
#             executor.submit(compute_params, input_prefix, param)
