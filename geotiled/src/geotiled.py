'''
GEOtiled Refactored Library v1.0.0
GCLab 2024

Compiled by Jay Ashworth (@washwor1) and Gabriel Laboy (@glaboy-vol)

Derived from original work by: Camila Roa (@CamilaR20), Eric Vaughan (@VaughanEric), Andrew Mueller (@Andym1098), Sam Baumann (@sam-baumann), David Huang (@dhuang0212), and Ben Klein (@robobenklein)

Learn more about GEOtiled from the paper: https://dl.acm.org/doi/pdf/10.1145/3588195.3595941
'''

from osgeo import osr, ogr, gdal
from datetime import datetime
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
import time
import os

# To install in Ubuntu: 1) sudo apt-get install grass grass-doc 2) pip install grass-session
# Need to do some work to get extensions with scripts
# https://github.com/OSGeo/grass-addons (repo with addon code and installation)
# https://grass.osgeo.org/grass83/manuals/g.extension.html (how to use g.extension)
# https://grass.osgeo.org/grass83/manuals/addons/r.valley.bottom.html (about valley depth script)
# https://grasswiki.osgeo.org/wiki/GRASS-QGIS_relevant_module_list (default script list - look at scripts starting with 'r.')
# Look in this directory: /usr/lib/grass78/etc/python/grass/script
from grass_session import Session
import grass.script as gscript
import tempfile

# CONSTANTS
TEXT_FILE_EXTENSION = ".txt"
SAGA_FILE_EXTENSION = ".sdat"
GEOTIFF_FILE_EXTENSION = ".tif"
VRT_DEFAULT_FILE_NAME = "merged.vrt"
SHAPEFILE_FOLDER_NAME = "shapefiles"


SHAPE_FILE_FOLDER_NAME = "shape_files"
DEM_FILE_FOLDER_NAME = "dem_tiles"

# USGS dataset codes used for fetch_dem function
USGS_DATASET_CODES = {"30m":"National Elevation Dataset (NED) 1 arc-second Current",
                      "10m":"National Elevation Dataset (NED) 1/3 arc-second Current"}

# GRASS modules needed for the parameters collected from https://grasswiki.osgeo.org/wiki/Terrain_analysis
# Look at r.stream.* for channel network stuff https://grasswiki.osgeo.org/wiki/Hydrological_Sciences
# For relative slope position: https://grass.osgeo.org/grass83/manuals/addons/r.slope.direction.html
TERRAIN_PARAMETER_CODES = {"SLP":"slope", 
                           "ASP":"aspect", 
                           "HLD":"hillshade", 
                           "CNL":"channel_network_base_level",
                           "CND":"channel_network_distance",
                           "CD":"closed_depressions",
                           "CI":"convergence_index",
                           "LSF":"ls_factor",
                           "PLC":"plan_curvature",
                           "PFC":"profile_curvature", 
                           "RSP":"relative_slope_position", 
                           "TCA":"total_catchment_area", 
                           "TWI":"topographic_wetness_index"}
                           "VD":"valley_depth"}

GRASS_PARAMETER_MODULES = {"slope":"r.slope.aspect", 
                           "aspect":"r.slope.aspect", 
                           "hillshade":"r.shade", # unsure 
                           "channel_network_base_level":"r.stream.channel", # unsure
                           "channel_network_distance":"r.stream.distance", # unsure
                           "closed_depressions":"r.fill.dir", # has two outputs https://grass.osgeo.org/grass83/manuals/r.fill.dir.html
                           "convergence_index":"r.convergence",
                           "ls_factor":"r.watershed",
                           "plan_curvature":"r.slope.aspect",
                           "profile_curvature":"r.slope.aspect", 
                           "relative_slope_position":"r.slope.direction", #unsure
                           "total_catchment_area":"r.catchment", # probably
                           "topographic_wetness_index":"r.topidx",
                           "valley_depth":"r.valley.bottom"}

SAGA_PARAMETER_FLAGS = {"elevation" : "ELEVATION",
                        "hillshade" : "SHADE",
                        "slope" : "SLOPE",
                        "aspect" : "ASPECT",
                        "plan_curvature" : "HCURV",
                        "profile_curvature" : "VCURV",
                        "convergence_index" : "CONVERGENCE",
                        "closed_depressions" : "SINKS",
                        "total_catchment_area" : "FLOW",
                        "topographic_wetness_index" : "WETNESS",
                        "ls_factor" : "LSFACTOR",
                        "channel_network_base_level" : "CHNL_BASE",
                        "channel_network_distance" : "CHNL_DIST",
                        "valley_depth" : "VALL_DEPTH",
                        "relative_slope_position" : "RSP"}

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

# Used to silence a deprecation warning. 
gdal.UseExceptions()

def __get_codes(code_type):
    """
    Reads code data from selected file and returns data as JSON dictionary.

    Parameters
    ----------
    code_type : str
        Specifies codes for which data are desired.
    """

    # Build path to text file containing region codes
    path_suffix = "codes/" + code_type + "_codes.txt"
    codes_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), path_suffix)

    # Read in JSON formatted codes
    data = None
    with open(codes_path) as f: 
        data = f.read() 
    codes = json.loads(data) 

    return codes

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def print_data_codes():
    """
    Print contents of 'data_codes.txt'.

    Outputs codes and their correlating terrain paramter from the 'data_codes.txt' file.
    These codes are meant to inform the user what elevation data is available for download off the USGS webpage.
    """

    # Read in JSON formatted codes
    data_codes = __get_codes("data")

    # Print all codes and their associated parameter
    print("The data codes and their associated parameter are:")
    for key in data_codes:
        print(key, ":", data_codes[key])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def print_parameter_codes(param=None):
    """
    Print contents of 'parameter_codes.txt'.

    Outputs codes and their correlating terrain paramter from the 'parameter_codes.txt' file.
    These codes are meant to inform the user what parameters are available for computation.

    Parameters
    ----------
    param : str, optional
        A terrain parameter to return the code for (default is None).
    """

    # Read in JSON formatted codes
    param_codes = __get_codes("parameter")

    # Output requested info
    if param is not None:
        param = param.lower()
        if param not in param_codes.values():
            # If the parameter doesn't exist, inform the user
            print("The entered parameter was not found.")
        else:
            # Print parameter and correlating code
            code_list = list(param_codes.keys())
            param_pos = list(param_codes.values()).index(param)
            print("The code for", param, "is", code_list[param_pos])
    else:
        # Print all codes and their associated parameter
        print("The parameter codes and their associated parameter are:")
        for key in param_codes:
            print(key, ":", param_codes[key])

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __extract_region_from_link(link):
    """
    Extract name of region from download link.

    The function pulls the name of the region from USGS download links for shapefiles.

    Parameters
    ----------
    link : str
        USGS specific link that specifies URL to download a shapefile.

    Returns
    -------
    str
        Region name extracted from the link.
    """
    
    re_link = "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_(.*)_State_Shape.zip"
    region = re.search(re_link, link).group(1).replace("_", " ") # Extract region from link
    return region

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def print_region_codes(region=None):
    """
    Print contents of 'region_codes.txt'.

    Outputs codes and their correlating region from the 'region_codes.txt' file.
    These code are meant to inform the user what shapefiles are available for download.

    Parameters
    ----------
    region : str, optional
        A region to return the code for (default is None). Not case-sensitive.
    """

    # Read in JSON formatted codes
    region_codes = __get_codes("region")

    # Output requested info
    if region is not None:
        region = region.lower()

        # Extract all regions and put them in a temp dictionary
        extracted_regions = {}
        for reg in region_codes.values(): 
            extracted_regions.update({__extract_region_from_link(reg).lower():reg})
            
        if region not in extracted_regions.keys():
            # If the code doesn't exist, inform the user
            print("The entered code was not found.")
        else:
            # Extract region from link and print
            code_list = list(region_codes.keys())
            reg_pos = list(region_codes.values()).index(extracted_regions[region])
            print("The code for", region, "is", code_list[reg_pos])
    else:
        # Print all codes and their associated region
        print("The region codes and their associated region are:")
        for key in region_codes:
            region = __extract_region_from_link(region_codes[key])
            print(key, ":", region)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __bash(argv):
    """
    Executes a command in bash.

    This function acts as a wrapper to execute bash commands using the Python subprocess Popen method. 
    Commands are executed synchronously, stdout and stderr are captured, and errors can be raised.

    Parameters
    ----------
    argv : list
        List of arguments for a bash command. They should be in the order that you would arrange them in the command line (e.g., ["ls", "-lh", "~/"]).

    Raises
    ------
    RuntimeError
        Popen returns with an error if the passed bash function returns an error.
    """

    arg_seq = [str(arg) for arg in argv] # Convert all arguments in list into a string
    proc = subprocess.Popen(arg_seq, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() # Synchronize
    stdout, stderr = proc.communicate()

    # Print error message if execution returned error
    if proc.returncode != 0:
        raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
            ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def set_working_directory(path):
    """
    Changes the working directory.

    Allows the user to change the working directory where data generated by GEOtiled will be stored.
    It will also create the directory for the user if it doesn't already exist, but the base 
    directory path must exist.

    Parameters
    ----------
    path : str
        Working directory to store data in.
    """

    # Ensure base directory path exists
    if os.path.exists(os.path.dirname(path)):
        Path(path).mkdir(parents=True, exist_ok=True)
        os.chdir(path) 
    else:
        print("Base path", os.path.dirname(path), "doesn't exist. Cannot set working directory.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
    
def __get_file_size(url):
    """
    Retrieve the size of a file at a given URL in bytes.

    This function sends a HEAD request to the provided URL and reads the 'Content-Length' header to determine the size of the file. 
    It's primarily designed to support the `download_files` function to calculate download sizes beforehand.

    Parameters
    ----------
    url : str
        The download URL to determine the file size from.

    Returns
    -------
    int
        Size of the file at the specified URL in bytes. Returns 0 if the size cannot be determined.
    """
    
    try:
        response = requests.head(url)
        return int(response.headers.get("Content-Length", 0))
    except:
        return 0

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __download_file(url, output_folder, pbar):
    """
    Download a single file from a URL and store it in a specified folder.

    This is a utility function that facilitates the downloading of files, especially within iterative download operations.
    If the file already exists, skip the download.

    Parameters
    ----------
    url : str
        The download URL of the file.
    output_folder : str
        The path where the downloaded file will be stored.
    pbar : tqdm object
        Reference to the tqdm progress bar, typically used in a parent function to indicate download progress.

    Returns
    -------
    int
        The number of bytes downloaded.
    """
    
    output_path = os.path.join(output_folder, url.split("/")[-1])
    if os.path.exists(output_path):
        return 0

    response = requests.get(url, stream=True)
    downloaded_size = 0
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
            downloaded_size += len(chunk)
            pbar.update(len(chunk))
            
    return downloaded_size

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_files(download_list, download_folder):
    """
    Download file(s) from provided URL(s) and store them in a desired folder.

    This function allows for the simultaneous downloading of files using threading
    and showcases download progress via a tqdm progress bar.

    Parameters
    ----------
    download_list : str, List[str]
        The name of a text file. This file should contain download URLs separated by newlines.
        A list of download URLs.
    download_folder : str
        Name of the folder where the downloaded files will be stored.
    """
    
    # Create folder to store downloaded files in
    output_folder = os.path.join(os.getcwd(), download_folder)
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Determine if URLs are raw list or from a text file
    if type(download_list) is not list:
        download_list = os.path.join(os.getcwd(), download_list + TEXT_FILE_EXTENSION) # Update list to be full path to text file
        with open(download_list, "r", encoding="utf8") as dsvfile:
            urls = [url.strip().replace("'$", "") for url in dsvfile.readlines()]
    else:
        urls = download_list
    
    # Compute total size of files to download
    total_size = sum(__get_file_size(url.strip()) for url in urls)

    # Create progress bar to track completion of downloads
    with tqdm(total=total_size, unit="B", unit_scale=True, ncols=100, desc="Downloading", colour="green") as pbar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            futures = [executor.submit(__download_file, url, output_folder, pbar) for url in urls]
            for future in concurrent.futures.as_completed(futures):
                size = future.result()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __download_shapefiles(codes):
    """
    Download shapefile(s) from a list of specified codes.
    
    This function downloads specified shapefile(s) off the USGS webpage and stores them in a 'shapefiles' folder located in the data directory.
    All codes passed should be valid codes from the 'region_codes.txt' file.
    It will skip downloading already existing shapefiles.

    Parameters
    ----------
    codes : str, List[str]
        Comma-separated string of valid codes.
        List of valid codes.
    """
    
    # Turn codes into list if not one
    if not isinstance(codes, list): 
        codes = list(codes.split(","))
    
    # Ensure codes passed are valid codes
    region_codes = __get_codes("region")
    for code in codes:
        if code not in list(region_codes.keys()):
            print(code, "is not a valid region code. Terminating execution.")
            return
    
    # Create path to shape file folder
    shape_path = os.path.join(os.getcwd(), SHAPEFILE_FOLDER_NAME)

    # Iterate through all passed code
    download_links = []
    download_codes = []
    for code in codes:
        # Create path to shape file
        shapefile_path = os.path.join(shape_path, code, code + ".shp")

        # Check if doesn't exist, add download URL to list
        if not os.path.isfile(shapefile_path):
            # Get download URL based off region code
            download_links.append(region_codes[code])
            download_codes.append(code)

    # Do the downloads if the list contains URLs
    if len(download_links) > 0:
        download_files(download_links, SHAPEFILE_FOLDER_NAME)

        # Unzip and rename files
        for iter, link in enumerate(download_links):
            # Isolate zip file name
            zip_name = os.path.basename(link)

            # Update path to zip file
            zip_path = os.path.join(shape_path, zip_name)

            # Unzip and get correct region files
            original_file_name = "GU_StateOrTerritory"
            original_file_path = "Shape/" + original_file_name
            with zipfile.ZipFile(zip_path, "r") as file:
                file.extract(original_file_path+".dbf", path=shape_path)
                file.extract(original_file_path+".prj", path=shape_path)
                file.extract(original_file_path+".shp", path=shape_path)
                file.extract(original_file_path+".shx", path=shape_path)
            
            file.close()
            
            # Rename folder and files
            dwnld_code = download_codes[iter]
            new_directory = os.path.join(shape_path, dwnld_code)
            new_dir_with_orig = os.path.join(new_directory, original_file_name)
            new_dir_with_code = os.path.join(new_directory, dwnld_code)
            os.rename(os.path.join(shape_path, "Shape"), new_directory)
            os.rename(new_dir_with_orig+".dbf", new_dir_with_code+".dbf")
            os.rename(new_dir_with_orig+".prj", new_dir_with_code+".prj")
            os.rename(new_dir_with_orig+".shp", new_dir_with_code+".shp")
            os.rename(new_dir_with_orig+".shx", new_dir_with_code+".shx")

            # Delete original zip file
            os.remove(zip_path)
            
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __get_extents(shapefile):
    """
    Returns the bounding box (extents) of a shapefile.
    
    This function extracts the extents of a shapefile. The extent is the upper left and lower right coordinates of the file.
    If the shapefile passed in isn't already download it, this function will do it automatically.

    Parameters
    ----------
    shapefile : str
        Code of shapefile to get extents for.
    
    Returns
    -------
    Tuple(Tuple(float))
        Returns two tuples - the first being the upper left (x,y) coordinate, and the second being the lower right (x,y) coordinate
    """
    
    # Build file path to shapefile, and if it doesn't exist go ahead and download it
    shape_path = os.path.join(os.getcwd(), SHAPEFILE_FOLDER_NAME, shapefile, shapefile + ".shp")
    if not os.path.exists(shape_path):
        __download_shapefiles(shapefile)
    
    ds = ogr.Open(shape_path)
    layer = ds.GetLayer()
    ext = layer.GetExtent()
    upper_left = (ext[0], ext[3])
    lower_right = (ext[1], ext[2])

    return upper_left, lower_right 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def fetch_dem(shapefile=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset="30m", prod_format="GeoTIFF", txt_file="download_urls", save_to_txt=True, download_folder="dem_tiles", download=False):
    """
    Queries USGS API for DEM data (URLs) given specified parameters.

    This function targets USGS National Map API to fetch DEM data URLs using specified parameters and can either 
    save the list of URLs to a text file or download from the list of URLs immediately. 
    GEOtiled currently only supports computation on GeoTIFF files, so it is not recommended to change the `prod_format` variable.

    Parameters
    ----------
    shapefile : str
        Code of shapefile with which a bounding box will be generated (default is None). Overrides the 'bbox' parameter if set.
    bbox : dict
        Bounding box coordinates to query for DEM data (default is {"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}).
    dataset : str
        Code for DEM data to download (default is 30m).
    prod_format : str
        File type of DEM data to download (default is GeoTIFF).
    txt_file : str
        Name of text file to save URLs to (default is download_urls).
    save_to_txt : bool
        Allow DEM URLs to be saved to a text file (default is True).
    download_folder : str
        Name of the download folder (default is dem_tiles).
    download : bool
        Allow DEM URLs retrieved to be immediately downloaded (default is False).

    Raises
    ------
    Exception
        If an invalid response while getting URLs off the USGS webpage, the function will error.
    """
    
    # Get coordinate extents if a shape file was specified
    if shapefile is not None:
        print("Reading in shape file...")

        # Get extents of shape file
        coords = __get_extents(shapefile)
        bbox["xmin"] = coords[0][0]
        bbox["ymax"] = coords[0][1]
        bbox["xmax"] = coords[1][0]
        bbox["ymin"] = coords[1][1]

    # Construct the query parameters
    print("Setting boundary extents...")
    data_codes = __get_codes("data")
    params = {
        "bbox": f"{bbox["xmin"]},{bbox["ymin"]},{bbox["xmax"]},{bbox["ymax"]}",
        "datasets": data_codes[dataset],
        "prodFormats": prod_format
    }

    # Make a GET request
    print("Requesting data from USGS...")
    base_url = "https://tnmaccess.nationalmap.gov/api/v1/products"
    response = requests.get(base_url, params=params)

    # Check for a successful request
    if response.status_code != 200:
        raise Exception(
            f"Failed to fetch data. Status code: {response.status_code}")

    # Convert JSON response to Python dict
    data = response.json()

    # Extract download URLs
    download_urls = [item["downloadURL"] for item in data["items"]]

    # Save URLs to text file 
    if save_to_txt is True:
        print("Saving URLs to text file...")
        txt_Path = os.path.join(os.getcwd(), txt_file + TEXT_FILE_EXTENSION) # Build full path to text file
        with open(txt_Path, "w") as file:
            for url in download_urls:
                file.write(f"{url}\n")

    # Download the files from the URLs
    if download is True:
        download_files(download_urls, download_folder)

    print("Fetch process complete.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_mosaic(input_folder, output_file, description, cleanup=False):
    """
    Builds a mosaic out of multiple GeoTIFF files.
    
    This function creates a mosaic from a list of GeoTIFF files utilizing the buildVRT() function from the GDAL library.

    Parameters
    ----------
    input_folder : str
        Name of the folder in data directory where files to mosaic together are located.
    output_file : str
        Name of mosaic file produced.
    description : str
        Description to add to output raster band of mosaic file.
    cleanup : bool, optional
        Determines if files from `input_folder` should be deleted after mosaic is complete (default is False).
    """
    
    # Create paths to files needed for computation
    vrt_path = os.path.join(os.getcwd(), VRT_DEFAULT_FILE_NAME)
    mosaic_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    input_folder = os.path.join(os.getcwd(), input_folder)
    
    # Check for valid folder and put input files into list
    print("Getting input files...")
    if not Path(input_folder).exists():
        print("The folder", input_folder, "does not exist. Terminating execution.")
        return
    input_files = glob.glob(input_folder + "/*" + GEOTIFF_FILE_EXTENSION)

    # Build VRT (mosaic)
    print("Constructing VRT...")
    vrt = gdal.BuildVRT(vrt_path, input_files)
    translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"])
    gdal.Translate(mosaic_path, vrt, options=translate_options)
    vrt = None  # close file

    # Update band description with name of terrain parameter
    print("Updating band description...")
    dataset = gdal.Open(mosaic_path)
    band = dataset.GetRasterBand(1)
    band.SetDescription(description)
    dataset = None  # close file

    # Delete intermediary tiles used to build mosaic
    if cleanup is True:
        print("Cleaning intermediary files...")
        shutil.rmtree(input_folder)
    os.remove(vrt_path)

    print("Mosaic process complete.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def reproject(input_file, output_file, projection, cleanup=False):
    """
    Reprojects a GeoTIFF file to a specified projection.
    
    This function reprojects a specified GeoTIFF file to a new projection changing the coordinate system representation, 
    and the result is saved to a new file. Multithreading is utilized to improve performance.

    Parameters
    ----------
    input_file : str
        Name of GeoTIFF file in data directory to reproject.
    output_file : str
        Name of new reprojected GeoTIFF file to store in data directory.
    projection : str
        Projection to use for reprojection. Can be a EPSG code (e.g. EPSG:4326) or the path to a WKT file.
    cleanup : bool, optional
        Determines if `input_file` should be deleted after reprojection is complete (default is False).
    """
    
    # Set full file paths for input and output
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)  

    # Ensure file to reproject exists
    if not os.path.isfile(input_path):
        print(os.path.basename(input_path), "does not exist. Terminating execution.")
        return
    
    print("Reprojection of", os.path.basename(input_path), "has begun.")
    
    # Set warp options for reprojection and warp
    warp_options = gdal.WarpOptions(dstSRS=projection, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"],
                                    multithread=True, warpOptions=["NUM_THREADS=ALL_CPUS"])
    warp = gdal.Warp(output_path, input_path, options=warp_options)
    warp = None  # Close file

    # Delete files used for reprojection
    if cleanup is True:
        print("Cleaning intermediary files...")
        os.remove(input_path)
        os.remove(input_path + ".aux.xml")

    print("Reprojection process complete.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __crop_pixels(input_file, output_file, window):
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
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)  
    
    # Window to crop by [left_x, top_y, width, height]
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

    # Perform translation
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_into_tiles(input_file, output_folder, num_tiles, buffer=10):
    """
    Splits a GeoTIFF file into smaller, equally-sized tiles.
    
    This function divides a GeoTIFF file into a specified number of tiles with added buffer regions to assist with rapid 
    computation of parameters by GEOtiled.

    Parameters 
    ----------
    input_file : str
        Name of the GeoTIFF file in the data directory to crop.
    output_folder : str
        Name of the folder in the data directory to store the cropped tiles.
    n_tiles : int 
        Number of total tiles to produce. Should be a perfect square number.
    buffer : int
        Specifies the buffer size - overlapping pixels that is included in the borders between two tiles (default is 10).
    """
    
    # Update path to input file
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)

    # Ensure file to crop exists
    if not os.path.isfile(input_path):
        print(os.path.basename(input_path), "does not exist. Terminating execution.")
        return
    
    # Update path to out folder and create it if not done so
    output_path = os.path.join(os.getcwd(), output_folder)
    Path(output_path).mkdir(parents=True, exist_ok=True)
    
    # Square root number of tiles to help get even number of rows and columns
    num_tiles = math.sqrt(num_tiles)

    # Split rows and columns of original file into even number of pixels for total number of tiles specified
    ds = gdal.Open(input_path, 0)
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

            tile_file = output_path + "/tile_" + "{0:04d}".format(tile_count) + GEOTIFF_FILE_EXTENSION
            win = [j, i, ncols, nrows]

            # Upper left corner
            win[0] = max(0, win[0] - buffer)
            win[1] = max(0, win[1] - buffer)

            w = win[2] + 2 * buffer
            win[2] = w if win[0] + w < cols else cols - win[0]

            h = win[3] + 2 * buffer
            win[3] = h if win[1] + h < rows else rows - win[1]  

            __crop_pixels(input_file, os.path.join(output_folder, os.path.basename(tile_file).replace(GEOTIFF_FILE_EXTENSION, "")), win)
            print(os.path.basename(tile_file), "cropped.")
            tile_count += 1 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __compute_params(input_path, param_list):
    """
    Computes terrain parameters for a given elevation GeoTIFF file.
    
    This function utilizes GDAL and GRASS libraries to compute parameters like slope, aspect, and more for a GeoTIFF file containing elevation data.
    The full parameter list can be found in the 'parameter_codes.txt' file.

    Parameters 
    ----------
    input_path : str
        Path to a GeoTIFF elevation file to compute parameters with.
    param_list : List[str]
        List of valid parameters to compute. The prefix for folder names should be the last item in the list.
    """
    
    # Get the name of the elevation file
    input_file = os.path.basename(input_path)

    # Compute all terrain parameters from passed in list
    for param in param_list:
        if param != param_list[-1]:
            # Set directory where computed tile will be stored
            param_path = os.path.join(os.getcwd(), param_list[-1] + param + "_tiles")
            Path(param_path).mkdir(parents=True, exist_ok=True)
            output_path = os.path.join(param_path, input_file)
            
            # Set correct options and compute parameter
            if param in ["slope", "aspect", "hillshade"]:
                if param == "aspect":
                    dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
                    gdal.DEMProcessing(output_path, input_path, processing=param, options=dem_options)
                else:
                    dem_options = gdal.DEMProcessingOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
                    gdal.DEMProcessing(output_path, input_path, processing=param, options=dem_options)
            else:
                # Define where to process the data in the temporary grass-session
                tmpdir = tempfile.TemporaryDirectory()
    
                # Create GRASS session
                s = Session()
                s.open(gisdb=tmpdir.name, location='PERMANENT', create_opts=input_path)
                creation_options = 'BIGTIFF=YES,COMPRESS=LZW,TILED=YES' # For GeoTIFF files
    
                # https://grasswiki.osgeo.org/wiki/GRASS_Python_Scripting_Library <- about the GRASS run_command function
                
                # Load raster into GRASS without loading it into memory (else use r.import or r.in.gdal)
                gscript.run_command('r.external', input=input_path, output='elevation', overwrite=True, quiet=True)
                
                # Set output folder for computed parameters
                gscript.run_command('r.external.out', directory=param_path, format="GTiff", option=creation_options, quiet=True)
    
                # Compute parameter
                if param == 'topographic_wetness_index':
                    gscript.run_command('r.topidx', input='elevation', output=input_file, overwrite=True, quiet=True)
                elif param == 'plan_curvature':
                    gscript.run_command('r.slope.aspect', elevation='elevation', tcurvature=input_file, flags='e', overwrite=True, quiet=True)
                elif param == 'profile_curvature':
                    gscript.run_command('r.slope.aspect', elevation='elevation', pcurvature=input_file, flags='e', overwrite=True, quiet=True)
                elif param == 'convergence_index':
                    gscript.run_command('g.extension', extension='r.convergence')
                    gscript.run_command('r.convergence', input='elevation', output=input_file, overwrite=True, quiet=True) #addon
                elif param == 'valley_depth':
                    gscript.run_command('g.extension', extension='r.valley.depth')
                    gscript.run_command('r.valley.bottom', elevation='elevation', mrvbf=input_file, overwrite=True, quiet=True) #addon
                elif param == 'ls_factor':
                    gscript.run_command('r.watershed', elevation='elevation', threshold=1, length_slope=input_file, overwrite=True, quiet=True)
                
                # Cleanup
                tmpdir.cleanup()
                s.close()
        
            # Update band description (and nodata value for GRASS-computed params)
            dataset = gdal.Open(output_path)
            band = dataset.GetRasterBand(1)
            band.SetDescription(param)
            if param not in ["slope", "aspect", "hillshade"]:
                band.SetNoDataValue(-9999)
            dataset = None

    print("Computation of parameters for", input_file.replace(GEOTIFF_FILE_EXTENSION,""), "completed.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __tmux_session_check():
    """
    Determine if there are any active tmux sessions.

    Uses Python subprocess to send out a 'tmux list-sessions' command to determine if there are active tmux sessions.
    Determination is done by evaluating the stderr for the message 'no server running'.

    Returns
    -------
    bool
        Specifies if there are active tmux sessions or not.
    """

    # Run the Python subprocess
    proc = subprocess.Popen(["tmux", "list-sessions"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() # Synchronize
    stdout, stderr = proc.communicate()

    # Evaluate the stderr
    if "no server running" in str(stderr):
        return True
    else:
        return False

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __geotiff_to_sdat(input_path, output_path): 
    """
    Convert a GeoTIFF file to an SDAT file.

    Uses GDAL to convert a GDAL supported GeoTIFF file to SAGA supported SDAT, SGRD, and PRJ files.

    Parameters
    ----------
    input_path : str
        Path to GeoTIFF file to convert.
    output_path : str
        Path to store created file in.
    """
    
    # Do the conversion
    translate_options = gdal.TranslateOptions(format="SAGA")
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __sdat_to_geotiff(input_file, output_file): 
    """
    Convert an SDAT file to a GeoTIFF file.

    Uses GDAL to convert a SAGA supported SDAT file to a GDAL supported GeoTIFF file.

    Parameters
    ----------
    input_path : str
        Path to GeoTIFF file to convert.
    output_path : str
        Path to store created file in.
    """

    # Do the conversion
    translate_options = gdal.TranslateOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
    gdal.Translate(output_file, input_file, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __compute_params_saga(input_files, param_list):
    """
    Compute terrain parameters using SAGA GIS.

    Computes a list of specified terrain parameters using the SAGA 'saga_cmd ta_compound 0' command.
    The command is run in a tmux session, and the function completes once the tmux session is terminated.
    The function handles conversion of the input elevation files from GeoTIFF to SAGA supported SDAT, SGRD, and PROJ.
    Note that SAGA computes parameters substantially slower than GDAL or GRASS.

    Parameters
    ----------
    input_files : List[str]
        List of elevation file input paths to compute parameters for.
    param_list : List[str]
        List of valid terrain parameters to compute for.
    """

    print('Computing parameters using SAGA...')
    
    # Compute for each input file
    for file_path in input_files:
        # Get some path information about the elevation file
        elevation_file = os.path.basename(file_path)
        elevation_path = os.path.dirname(file_path)
        elevation_folder = elevation_path[elevation_path.rindex('/')+1:]

        # Convert names to store files in SAGA correlated folders
        saga_elevation_file = elevation_file.replace(GEOTIFF_FILE_EXTENSION, SAGA_FILE_EXTENSION)
        saga_elevation_path = elevation_path.replace(elevation_folder, "saga_" + elevation_folder)
        Path(saga_elevation_path).mkdir(parents=True, exist_ok=True)
        saga_file_path = os.path.join(saga_elevation_path, saga_elevation_file)

        # Convert GeoTIFF file to a SDAT file
        print("Converting", elevation_file, "to SDAT...")
        __geotiff_to_sdat(file_path, saga_file_path)
        
        # Initial command setup
        print("Building command...")
        cmd = ["tmux", "new-session", "-d", "-s", "terrainParamsSession", "saga_cmd", "ta_compound", "0", "-ELEVATION", saga_file_path.replace(SAGA_FILE_EXTENSION, ".sgrd")]
        
        # Add parameters to the command
        for param in param_list:
            # Create parameter folder if it doesnt already exist
            saga_elevation_folder = saga_elevation_path[saga_elevation_path.rindex('/')+1:]
            param_path = saga_elevation_path.replace(saga_elevation_folder, "saga_" + param + "_tiles")
            Path(param_path).mkdir(parents=True, exist_ok=True)
            param_file = elevation_file.replace(GEOTIFF_FILE_EXTENSION, ".sgrd")
            
            cmd.append("-" + SAGA_PARAMETER_FLAGS[param])
            cmd.append(os.path.join(param_path, param_file))

        # Run the command and wait for it to finish
        print("Computing parameters...")
        __bash(cmd)
        time.sleep(1)
        while not __tmux_session_check():
            time.sleep(3)

        print("Computation of parameters for", elevation_file, "complete.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_geotiled(input_folder, param_list, num_procs, output_folder_prefix='', use_saga=False, cleanup=False):
    """
    Configures the multiprocessing pool for GEOtiled to begin computing terrain parameters.
    
    This function utilizes the multiprocessing library to allow for computation of parameters on different elevation GeoTIFF files at the same time.
    Parameters that can be computed are specified in the 'parameter_codes.txt', or the 'all' keyword can be passed to compute all parameters.
    It is better to keep `num_procs` as a low value for systems with low amounts of RAM.
    
    Parameters
    ----------
    input_folder : str
        Name of the folder in the data directory containing DEM elevation files.
    param_list : List[str]
        List containing codes for terrain parameters to compute. The 'all' keyword will compute all params.
    num_procs : int
        Integer specifying the number of python instances to use for multiprocessing.
    output_folder_prefix : str
        Prefix to attach to all output folders created for storing computed terrain paramters (default is '').
    use_sage : bool
        Determine if parameters should be computed with SAGA (default is False).
    cleanup : bool
        Determine if elevation files should be deleted after computation (default is False).
    """
    
    # Check to ensure input folder exists
    input_path = os.path.join(os.getcwd(), input_folder)
    if not Path(input_path).exists():
        print("The folder", input_path, "does not exist. Terminating execution.")
        return
    
    # Check if all params are to be computed
    param_codes = __get_codes("parameter")
    if "all" in param_list:
        # Ensure only the 'all' code entered
        if len(param_list) > 1:
            print("Can't use 'all' keyword with other parameters in list")
            return
        
        param_list.remove("all")
        for key in param_codes:
            param_list.append(key)

    # Ensure all codes entered are valid codes and put full parameter name in different list
    params = []
    for param in param_list:
        if param not in list(param_codes.keys()):
            print("Parameter code", param, "is not a valid code")
            return
        params.append(param_codes[param])
    
    # Get files from input folder to compute on 
    print('Getting input files...')
    input_files = sorted(glob.glob(input_path + '/*' + GEOTIFF_FILE_EXTENSION))

    # Evaluate if params computed with SAGA or not
    if not use_saga:
        # Append the output_folder_prefix value to the params list to ensure it gets passed into the processing pool
        if output_folder_prefix != '':
            output_folder_prefix = output_folder_prefix + '_'
        params.append(output_folder_prefix)

        # Configure a list with the input files and selected params to support computing with pool.starmap()
        # https://superfastpython.com/multiprocessing-pool-starmap/
        items = []
        for input_file in input_files:
            items.append((input_file, params))

        # Start the multiprocessing pool
        print('Starting computation of parameters...')
        pool = multiprocessing.Pool(processes=num_procs) 
        pool.starmap(__compute_params, items)
    else:
        __compute_params_saga(input_files, params)

    # Remove files used to compute params
    if cleanup is True:
        print('Cleaning files...')
        shutil.rmtree(input_path)

    print('GEOtiled computation done!')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_mosaic_filtered(input_folder, output_file, cleanup=False):
    '''
    Builds mosaic from multiple GeoTIFF files that were cropped with buffer regions.
    --------------------------------------------------------------------------------
    This function is similar to the build_mosaic function but handles mosaicking together GeoTIFF files that were split to includes buffer regions and averages points in the buffer regions together when merging.

    Required Parameters
    -------------------
    input_files : str
        String specifying name of folder in data directory where GeoTIFF files to mosaic together are located.
    output_file : str
        String specifying name of mosaicked file produced.

    Optional Parameters
    -------------------
    cleanup : bool
        Boolean specifying if tiles in input_folder used for computation should be deleted after computation is complete. Default is False.

    Outputs
    -------
    GeoTIFF File
        Produces a GeoTIFF file that is the mosaic of all GeoTIFF files from the specified input folder.

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
    input_path = os.path.join(os.getcwd(), input_folder)
    if not Path(input_path).exists():
        print('The folder ' + input_path + ' does not exist. Terminating execution.')
        return
    
    print('Mosaicking started for ' + output_file + GEOTIFF_FILE_EXTENSION + '...')
    
    # Get files from input folder to merge together
    input_files = glob.glob(input_path + '/*' + GEOTIFF_FILE_EXTENSION)

    # Build full path for VRT and output file
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    vrt_file = os.path.join(os.getcwd(), 'merged.vrt')

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
           'BIGTIFF=YES', '--config', 'GDAL_VRT_ENABLE_PYTHON', 'YES', vrt_file, output_path]
    __bash(cmd)

    # Remove intermediary files used to build mosaic
    if cleanup is True:
        print('Cleaning files...')
        shutil.rmtree(input_path)
    os.remove(vrt_file)

    print('Mosaic process completed.')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_to_valid_data(input_file, output_file, block_size=512):
    '''
    Crops a border region of NaN values from a GeoTIFF file.
    --------------------------------------------------------
    This function uses blocking to scan through a GeoTIFF file to determine the extent of valid data and crops excess Nan values at the borders of the spatial data.

    Required Parameters 
    --------------------
    input_file : str
        String specifying name of GeoTIFF file to crop.
    output_file : str
        String specifying name of GeoTIFF file to save cropped data to.

    Optional Parameters
    --------------------
    block_size : int
        Integer specifying the block size to use when computing extents. Default is 512.

    Outputs
    -------
    GeoTIFF File
        Output is a cropped GeoTIFF file with the path specified by 'output_file'

    Returns
    -------
    None
        The function does not return any value.

    Notes
    ------
    - Blocking helps minimize RAM usage, and should be adjusted as needed to help improve performance.
    '''
    # Update file names to full paths
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    
    src_ds = gdal.Open(input_path, gdal.GA_ReadOnly)
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
    
    out_ds = gdal.Translate(output_path, src_ds, projWin=[min_x, max_y, max_x, min_y], projWinSRS='EPSG:4326')
    
    out_ds = None
    src_ds = None

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_region(codes, input_file, output_file):
    '''
    Crops a raster file based off the region defined by a combination of shapefiles.
    --------------------------------------------------------------------------------
    This function uses a GDAL function to crop a raster file according to boundaries specified by a shapefile.

    Required Parameters
    -------------------
    codes : str list
        String list specifying shapefile codes that will outline the cropping region.
    input_file : str
        String specifying name of input raster file to crop.
    output_file : str
        String specifying name of cropped output raster file to save.

    Outputs
    -------
    GeoTIFF File
        A GeoTIFF File with data cropped to the bounds of the shapefile is created.

    Returns
    -------
        The function returns a list with all paths to shapefiles used.

    Error states
    ------------
    - Error will generate if bounds of shapefile exceed input raster's limits.
    '''
    # Update file names to full paths
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)

    # Download shape files if any are missing
    __download_shapefiles(codes)

    # Update path to shape files
    shape_paths = []
    for code in codes:
        shape_paths.append(os.path.join(os.getcwd(), SHAPE_FILE_FOLDER_NAME, code, code + '.shp'))

    # Read each shapefile, clean any invalid geometries, and union them
    gdfs = [gpd.read_file(shp_file).buffer(0) for shp_file in shape_paths]
    combined_geom = gdfs[0].unary_union
    for gdf in gdfs[1:]:
        combined_geom = combined_geom.union(gdf.unary_union)
    
    combined_gdf = gpd.GeoDataFrame(geometry=[combined_geom], crs=gdfs[0].crs)
        
    # Save the combined shapefile temporarily for cropping
    temp_combined_shp = os.path.join(os.getcwd(), "temp_combined.shp")
    combined_gdf.to_file(temp_combined_shp)

    # Do the cropping
    warp_options = gdal.WarpOptions(cutlineDSName=temp_combined_shp, cropToCutline=True, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
    warp = gdal.Warp(output_path, input_path, options=warp_options)
    warp = None

    # Remove the temporary combined shapefile
    os.remove(temp_combined_shp)
    os.remove(os.path.join(os.getcwd(), "temp_combined.cpg"))
    os.remove(os.path.join(os.getcwd(), "temp_combined.dbf"))
    os.remove(os.path.join(os.getcwd(), "temp_combined.prj"))
    os.remove(os.path.join(os.getcwd(), "temp_combined.shx"))

    return shape_paths

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def generate_img(tif, cmap='inferno', dpi=150, downsample=1, verbose=False, clean=False, title=None, 
                 nancolor='green', ztype="Z", zunit=None, xyunit=None, vmin=None, vmax=None, reproject_gcs=False, 
                 shp_files=None, crop_shp=False, bordercolor="black", borderlinewidth=1.5, saveDir = None):
    '''
    Plots a GeoTIFF image using matplotlib.
    --------------------------------------
    This function is meant to visualize GeoTIFF with optional parameters to ensure data visualization is made easily discernable and customizable.

    Required Parameters
    --------------------
    tif : str
        String specifying name of the GeoTIFF file in the data directory to plot.

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
    vmin : float
        Value of lower bound for coloring on plot. Will be whatever the min value of the data is if none is specified.
    vmax : float
        Value of upper bound for coloring on plot. Will be whatever the max value of the data is if none is specified.
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
    saveDir : str
        String specifying directory to save image to.

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
    - Alternative colormaps can be found in the matplotlib documentation.
    - Graph currently only tested for visualization with Jupyter Notebooks.
    - Using 'shp_files' without setting 'crop_shp' will allow you to plot the outline of the shapefile without actually cropping anything.
    '''
    # Initial setup
    tif_dir_changed = False

    # Update full path to tif file
    tif_path = os.path.join(os.getcwd(), tif + GEOTIFF_FILE_EXTENSION)

    # Ensure file to plot exists
    if not os.path.isfile(tif_path):
        print(os.path.basename(tif_path) + ' does not exist. Terminating execution.' )
        return
    
    # Reproject raster into geographic coordinate system if needed
    if reproject_gcs:
        reproject_path = os.path.join(os.getcwd(), 'vis' + GEOTIFF_FILE_EXTENSION)
        reproject(tif, 'vis', "EPSG:4326")
        if crop_shp is False:
            crop_path = os.path.join(os.getcwd(), 'vis_trim_crop' + GEOTIFF_FILE_EXTENSION)
            print("Cropping NaN values...")
            crop_to_valid_data('vis', 'vis_trim_crop')
            os.remove(reproject_path)
            tif_path = crop_path
        else:
            tif_path = reproject_path
        tif_dir_changed = True

    # Crop using shapefiles if needed
    shape_paths = []
    if crop_shp and shp_files:
        # Check if the list is not empty
        if not shp_files:
            print("Shapefile list is empty. Skipping shapefile cropping.")
        else:                
            print("Cropping with combined shapefiles...")
            region_crop_path = os.path.join(os.getcwd(), 'crop' + GEOTIFF_FILE_EXTENSION)
            shape_paths = crop_region(shp_files, os.path.basename(tif_path).replace(GEOTIFF_FILE_EXTENSION,''), 'crop')
            
            if tif_dir_changed:
                os.remove(tif_path)
            tif_path = region_crop_path
            tif_dir_changed = True

    print("Reading in tif for visualization...")
    dataset = gdal.Open(tif_path)
    band = dataset.GetRasterBand(1)

    geotransform = dataset.GetGeoTransform()
    spatial_ref = osr.SpatialReference(wkt=dataset.GetProjection())

    # Extract spatial information about raster
    proj_name = spatial_ref.GetAttrValue('PROJECTION')
    proj_name = proj_name if proj_name else "GCS, No Projection"
    data_unit = zunit or spatial_ref.GetLinearUnitsName()
    coord_unit = xyunit or spatial_ref.GetAngularUnitsName()
    z_type = ztype if band.GetDescription() == '' else band.GetDescription()

    # Output verbose message if specified
    if verbose:
        print(f"Geotransform:\n{geotransform}\n\nSpatial Reference:\n{spatial_ref}\n\nDocumentation on spatial reference format: https://docs.ogc.org/is/18-010r11/18-010r11.pdf\n")

    raster_array = gdal.Warp('', tif_path, format='MEM', 
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

    vmn = vmin
    vmx = vmax
    if vmin is None:
        vmn = np.nanmin(raster_array)
    if vmax is None:
        vmx = np.nanmax(raster_array)
        
    sm = ax.imshow(raster_array, cmap=cmap_instance, vmin=vmn, vmax=vmx,
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
        os.remove(tif_path)

    return raster_array

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_coord(input_file, output_file, upper_left, lower_right):
    '''
    Crops raster file to a specific region given specified coordinates.
    -------------------------------------------------------------------
    This function uses GDAL functions to crop a raster file based on a specified upper-left and lower-right coordinate.

    Required Parameters
    -------------------
    input_file : str
        String specifying name of the input raster file to crop.
    output_file : str
        String specifying name of cropped raster to save.
    upper_left : tuple of float
       Float tuple specifying upper-left (x,y) coordinates to crop raster from.
    lower_right : tuple of float
        Float tuple specifying lower-right (x,y) coordinates to crop raster from.

    Outputs
    -------
    GeoTIFF File
        Will generate the cropped GeoTIFF file at the location specified by 'output_file'.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - The `upper_left` and `lower_right` coordinates define the bounding box for cropping.
    - Coordinates must be in the same projection as the raster

    Error states
    ------------
    - An error will be thrown if coordinates specified fall outside the range of the input raster's bounds.
    '''
    # Build full paths
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    
    # Crop
    window = upper_left + lower_right
    translate_options = gdal.TranslateOptions(projWin=window, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])#,callback=gdal.TermProgress_nocb)
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def change_raster_format(input_file, output_file, raster_format):
    '''
    Convert format of a specified raster file.
    ------------------------------------------
    This function uses GDAL functions to convert the file type of specified raster data.

    Required Parameters
    -------------------
    input_file : str
        String specifying name of input GeoTIFF raster file.
    output_file : str
        String specifying name of output raster file.
    raster_format : str
        GDAL supported string specifying the raster format to convert the input file to.

    Outputs
    -------
    Raster File
        Generates a raster of the format specified by 'raster_format'.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - Supported formats can be found of the GDAL website: https://gdal.org/drivers/raster/index.html
    '''
    # Build full paths
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    
    # SAGA, GTiff
    if raster_format == 'GTiff':
        translate_options = gdal.TranslateOptions(format=raster_format, creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
    elif raster_format == 'NC4C':
        translate_options = gdal.TranslateOptions(format=raster_format, creationOptions=['COMPRESS=DEFLATE'])
    else:
        translate_options = gdal.TranslateOptions(format=raster_format)
    
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def extract_raster(csv_file, raster_file, band_names):
    '''
    Extracts raster values and stores them with correlating coordinates found in a CSV file.
    ----------------------------------------------------------------------------------------
    This function reads raster values from a file and stores the value with correlating x and y coordinates in a CSV file.

    Required Parameters
    -------------------
    csv_file : str
        String specifying name of CSV file for input.
    raster_file : str
        String specifying name to raster file to read raster values from.
    band_names : str
        String list specifying names of bands to extract values from in input raster.

    Outputs
    -------
    CSV File
        An updated CSV file with raster values correlating to x and y coordinates will be produced.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - CSV file must already be appropiately formatted with x and y coordinates.
    - Order of 'band_names' should correlate to order of band names in specified raster file.
    - If coordinates in CSV file are outside of bounds of raster, incorrect or no values will be extracted.

    Error States
    ------------
    - If the CSV file does not have 'x' and 'y' columns, a KeyError will occur.
    - If the specified coordinates in the CSV file are outside the bounds of the raster, incorrect or no values may be extracted.
    '''
    # Build full paths
    csv_path = os.path.join(os.getcwd(), csv_file + '.csv')
    raster_path = os.path.join(os.getcwd(), raster_file + GEOTIFF_FILE_EXTENSION)
    
    # Extract values from raster corresponding to
    df = pd.read_csv(csv_path)

    ds = gdal.Open(raster_path, 0)
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

    df.to_csv(csv_path, index=None)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_stack(input_folder, output_file):
    '''
    Stacks multiple GeoTIFF files into a single GeoTIFF with multiple bands.
    ------------------------------------------------------------------------
    This function takes multiple GeoTIFF files and combines them into one file represented by different bands. This is useful when multiple datasets need to be represented as one.

    Required Parameters
    -------------------
    input_folder : str
        String specifying name of folder to GeoTIFF files to be stacked together.
    output_file : str
        String specifying name of file to store stacked files to.

    Outputs
    -------
    GeoTIFF File
        Generates a single GeoTIFF file with multiple bands with the path it is stored in being establed by the variable 'output_file'.

    Returns
    -------
        The function does not return any value.

    Notes
    -----
    - The band order will be based off the 'input_files' list order
    '''
    # Build full paths
    input_path = os.path.join(os.getcwd(), input_folder)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)    

    print('Getting input files...')
    if not Path(input_path).exists():
        print('The folder ' + input_path + ' does not exist. Terminating execution.')
        return
    input_files = glob.glob(input_path + '/*' + GEOTIFF_FILE_EXTENSION)
    
    vrt_options = gdal.BuildVRTOptions(separate=True)
    vrt = gdal.BuildVRT("stack.vrt", input_files, options=vrt_options)
    translate_options = gdal.TranslateOptions(creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])#,callback=gdal.TermProgress_nocb)
    gdal.Translate(output_path, vrt, options=translate_options)
    vrt = None  # closes file
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def tif2csv(input_file, output_file='params', band_names=['elevation']):
    '''
    Converts a raster file into a CSV file.
    ---------------------------------------
    This function reads values from a GeoTIFF file and converts it into the CSV format. The CSV columns are the x coordinate, y coordinate, and raster value.

    Required Parameters
    -------------------
    input_file : str
        String specifying name of GeoTIFF file to be stored as a CSV.

    Optional Parameters
    -------------------
    output_file : str
        String specifying name of CSV file to save data to. Default is 'params'.
    band_names : list
        String list specifying names of GeoTIFF bands to pull data from. Default is ['elevation'].

    Outputs
    -------
    CSV File
        Generates CSV file with converted data from raster file.

    Returns
    -------
        The function does not return any value.

    Notes
    -----
    - NaN values indicate no data found at a particular coordinate
    - Only band names provided will be read into the CSV. If missing band names from GeoTIFF file, CSV may contain columns without headers or be missing data.
    '''
    # Build full paths
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + '.csv')

    ds = gdal.Open(input_path, 0)
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
    df.to_csv(output_path, index=None)

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
