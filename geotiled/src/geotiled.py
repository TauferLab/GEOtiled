"""
GEOtiled Library v0.2
GCLab 2024

Compiled by Jay Ashworth (@washwor1) and Gabriel Laboy (@glaboy-vol)

Derived from original work by: Camila Roa (@CamilaR20), Eric Vaughan (@VaughanEric), Andrew Mueller (@Andym1098), Sam Baumann (@sam-baumann), David Huang (@dhuang0212), and Ben Klein (@robobenklein)

Learn more about GEOtiled from the paper: https://dl.acm.org/doi/pdf/10.1145/3588195.3595941
"""

from osgeo import osr, ogr, gdal
from pathlib import Path
from tqdm import tqdm

import matplotlib.pyplot as plt
import saga_wrapper as sw
import geopandas as gpd
import pandas as pd
import numpy as np

import concurrent.futures
import multiprocessing
import subprocess
import warnings
import requests
import zipfile
import shutil
import glob
import json
import os
import re

# Enable exceptions for GDAL 
gdal.UseExceptions()

# Suppress specific warning related to loading shapefiles
warnings.filterwarnings("ignore", message=r"Measured \(M\) geometry types are not supported.*")

#################
### CONSTANTS ###
#################

SHAPEFILE_FOLDER_NAME = "shapefiles"
VRT_DEFAULT_FILE_NAME = "merged.vrt"

# List is computing individual parameters
COMPUTABLE_PARAMETERS = ["slope", "aspect", "hillshade", "plan_curvature", "profile_curvature", "convergence_index", "total_catchment_area", "specific_catchment_area", "topographic_wetness_index", "ls_factor", "channel_network", "drainage_basins", "channel_network_base_level", "channel_network_distance", "valley_depth", "relative_slope_position"]

# List if using 'all' parameters
SAGA_PARAMETER_LIST = ["slope", "aspect", "hillshade", "plan_curvature", "profile_curvature", "convergence_index", "closed_depressions", "total_catchment_area", "topographic_wetness_index", "ls_factor", "channel_network", "drainage_basins", "channel_network_base_level", "channel_network_distance", "valley_depth", "relative_slope_position"]

DATA_CODES = {"60m": "National Elevation Dataset (NED) Alaska 2 arc-second Current",
              "30m": "National Elevation Dataset (NED) 1 arc-second Current",
              "10m": "National Elevation Dataset (NED) 1/3 arc-second Current",
              "5m": "Alaska IFSAR 5 meter DEM",
              "3m": "National Elevation Dataset (NED) 1/9 arc-second",
              "1m": "Digital Elevation Model (DEM) 1 meter"}

DATA_FORMATS = {"60m": "GeoTIFF",
                "30m": "GeoTIFF",
                "10m": "GeoTIFF",
                "5m": "ArcGrid",
                "3m": "IMG",
                "1m": "GeoTIFF"}

REGION_CODES = {"AL": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Alabama_State_Shape.zip", 
                 "AK": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Alaska_State_Shape.zip", 
                 "AZ": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Arizona_State_Shape.zip", 
                 "AR": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Arkansas_State_Shape.zip", 
                 "CA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_California_State_Shape.zip", 
                 "CO": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Colorado_State_Shape.zip", 
                 "CT": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Connecticut_State_Shape.zip", 
                 "DE": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Delaware_State_Shape.zip", 
                 "DC": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_District_of_Columbia_State_Shape.zip",
                 "FL": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Florida_State_Shape.zip", 
                 "GA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Georgia_State_Shape.zip", 
                 "GU": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Guam_State_Shape.zip",
                 "HI": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Hawaii_State_Shape.zip", 
                 "ID": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Idaho_State_Shape.zip",
                 "IL": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Illinois_State_Shape.zip", 
                 "IN": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Indiana_State_Shape.zip", 
                 "IA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Iowa_State_Shape.zip", 
                 "KS": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Kansas_State_Shape.zip", 
                 "KY": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Kentucky_State_Shape.zip",
                 "LA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Louisiana_State_Shape.zip", 
                 "ME": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Maine_State_Shape.zip", 
                 "MD": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Maryland_State_Shape.zip", 
                 "MA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Massachusetts_State_Shape.zip", 
                 "MI": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Michigan_State_Shape.zip",
                 "MN": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Minnesota_State_Shape.zip", 
                 "MS": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Mississippi_State_Shape.zip", 
                 "MO": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Missouri_State_Shape.zip", 
                 "MT": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Montana_State_Shape.zip", 
                 "NE": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Nebraska_State_Shape.zip",
                 "NV": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Nevada_State_Shape.zip", 
                 "NH": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_Hampshire_State_Shape.zip", 
                 "NJ": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_Jersey_State_Shape.zip", 
                 "NM": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_Mexico_State_Shape.zip", 
                 "NY": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_New_York_State_Shape.zip",
                 "NC": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_North_Carolina_State_Shape.zip", 
                 "ND": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_North_Dakota_State_Shape.zip", 
                 "MP": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Commonwealth_of_the_Northern_Mariana_Islands_State_Shape.zip",
                 "OH": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Ohio_State_Shape.zip", 
                 "OK": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Oklahoma_State_Shape.zip",
                 "OR": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Oregon_State_Shape.zip", 
                 "PA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Pennsylvania_State_Shape.zip", 
                 "PR": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Puerto_Rico_State_Shape.zip", 
                 "RI": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Rhode_Island_State_Shape.zip", 
                 "SC": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_South_Carolina_State_Shape.zip",
                 "SD": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_South_Dakota_State_Shape.zip", 
                 "TN": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Tennessee_State_Shape.zip", 
                 "TX": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Texas_State_Shape.zip", 
                 "UT": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Utah_State_Shape.zip", 
                 "VT": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Vermont_State_Shape.zip",
                 "VA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Virginia_State_Shape.zip", 
                 "VI": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_United_States_Virgin_Islands_State_Shape.zip", 
                 "WA": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Washington_State_Shape.zip", 
                 "WV": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_West_Virginia_State_Shape.zip", 
                 "WI": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Wisconsin_State_Shape.zip",
                 "WY": "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_Wyoming_State_Shape.zip"}

###############################
### MISCELLANEOUS FUNCTIONS ###
###############################

def determine_if_path(path):
    """
    Determine if entered path is a folder or filename or full path.

    Determines if an entered path is a folder or filename or full path.
    If the entered path is a full, valid path, the path is returned.
    If the entered path is only a folder or filename, a path is built 
    with the specified name and current working directory and the 
    new full path is returned.

    Parameters
    ----------
    path : str
        A full path or folder/file name.
    """

    # Determine if the entered path is just a folder/filename 
    if Path(path).exists():
        if '/' in path: 
            return path
        elif (os.getcwd() in path): 
            return path
        else: 
            return os.path.join(os.getcwd(), path)
    else:
        new_path = os.path.join(os.getcwd(), path)
        return new_path

###############################
### CONFIGURATION FUNCTIONS ###
###############################

def print_computable_parameters():
    """
    Print all computable parameters.

    Outputs all computable parameters for GEOtiled from list in script.
    A special indicator '(G)' indicates it is computable with GDAL.
    """

    print("Computable Terrain Parameters (G == GDAL compatible):")
    for parameter in COMPUTABLE_PARAMETERS:
        if (parameter == "slope") or (parameter == "aspect") or (parameter == "hillshade"):
            print(parameter, "(G)")
        else:
            print(parameter)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def extract_region_from_link(link):
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

def print_region_codes():
    """
    Print contents of 'region_codes.txt'.

    Output codes and their correlating region from the REGION_CODES variable.
    These codes are meant to inform the user what shapefiles are available for download.
    """

    # Print all codes and their associated region
    print("Region codes and their associated region:")
    for key in REGION_CODES:
        region = extract_region_from_link(REGION_CODES[key])
        print(key, ":", region)

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

##################################
### INPUT VALIDATION FUNCTIONS ###
##################################

def validate_codes(input_codes, code_dictionary):
    for code in input_codes:
        if code not in list(code_dictionary.keys()):
            print(code, "is not a valid code. Terminating execution.")
            return -1
    return 0

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def validate_path_exists(input_path):
    if not Path(input_path).exists():
        print(input_path, "does not exist. Terminating execution.")
        return -1
    return 0

####################################
### FEATURE EXTRACTION FUNCTIONS ###
####################################

def get_extents(shapefile):
    """
    Returns the bounding box (extents) of a shapefile.
    
    This function extracts the extents of a shapefile. The extent is the upper left and lower right coordinates of the file.

    Parameters
    ----------
    shapefile : str
        Path to shapefile to get extents for.
    
    Returns
    -------
    Tuples(float)
        Returns two tuples - the first being the upper left (x,y) coordinate, and the second being the lower right (x,y) coordinate
    """

    # Get upper left and lower rights extents of shapefile and return
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer()
    ext = layer.GetExtent()
    upper_left = (ext[0], ext[3])
    lower_right = (ext[1], ext[2])

    return upper_left, lower_right 

###############################
### FILE DOWNLOAD FUNCTIONS ###
###############################

def get_file_size(url):
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

def download_file(url, output_folder, pbar):
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
        The name/path of a text file. This file should contain download URLs separated by newlines.
        A list of download URLs.
    download_folder : str
        Name/path of the folder where the downloaded files will be stored.
    """
    
    # Create folder to store downloaded files in
    output_folder = determine_if_path(download_folder)
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Determine if URLs are a list or from a text file
    if type(download_list) is not list:
        download_list = determine_if_path(download_list)
        with open(download_list, "r", encoding="utf8") as dsvfile:
            urls = [url.strip().replace("'$", "") for url in dsvfile.readlines()]
    else:
        urls = download_list
    
    # Compute total size of files to download
    total_size = sum(get_file_size(url.strip()) for url in urls)

    # Create progress bar to track completion of downloads
    with tqdm(total=total_size, unit="B", unit_scale=True, ncols=100, desc="Downloading", colour="green") as pbar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            futures = [executor.submit(download_file, url, output_folder, pbar) for url in urls]
            for future in concurrent.futures.as_completed(futures):
                size = future.result()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_shapefiles(codes):
    """
    Download specified shapefile(s) from a list of supplied codes.
    
    This function downloads specified shapefile(s) off the USGS webpage and stores them in a 'shapefiles' folder located in the working directory.
    All codes passed should be valid codes from the REGION_CODES variable.
    It will skip downloading already existing shapefiles.

    Parameters
    ----------
    codes : str, List[str]
        Comma-separated string of valid codes.
        List of valid codes.
    """
    
    # If codes are comma-separated, turn into a list
    if not isinstance(codes, list): 
        codes = list(codes.split(","))
    
    # Ensure codes passed are valid codes
    if validate_codes(codes, REGION_CODES) == -1: return
    
    # Create path to shapefile folder
    shapefile_folder_path = os.path.join(os.getcwd(), SHAPEFILE_FOLDER_NAME)

    # Iterate through all codes and download needed shapefiles
    download_links = []
    download_codes = []
    for code in codes:
        # Create path to shapefile
        shapefile_path = os.path.join(shapefile_folder_path, code, code + ".shp")

        # Check if shapefile not already downloaded, and add download URL to list
        if not os.path.isfile(shapefile_path):
            # Get download URL based off region code
            download_links.append(REGION_CODES[code])
            download_codes.append(code)

    # Begin download
    if len(download_links) > 0:
        download_files(download_links, SHAPEFILE_FOLDER_NAME)

        # Unzip and rename files
        for iter, link in enumerate(download_links):
            # Isolate zip file name
            zip_name = os.path.basename(link)

            # Update path to zip file
            zip_path = os.path.join(shapefile_folder_path, zip_name)

            # Unzip and get correct region files
            original_file_name = "GU_StateOrTerritory"
            original_file_path = "Shape/" + original_file_name
            with zipfile.ZipFile(zip_path, "r") as file:
                file.extract(original_file_path+".dbf", path=shapefile_folder_path)
                file.extract(original_file_path+".prj", path=shapefile_folder_path)
                file.extract(original_file_path+".shp", path=shapefile_folder_path)
                file.extract(original_file_path+".shx", path=shapefile_folder_path)
            
            file.close()
            
            # Rename folder and files
            dwnld_code = download_codes[iter]
            new_directory = os.path.join(shapefile_folder_path, dwnld_code)
            new_dir_with_orig = os.path.join(new_directory, original_file_name)
            new_dir_with_code = os.path.join(new_directory, dwnld_code)
            os.rename(os.path.join(shapefile_folder_path, "Shape"), new_directory)
            os.rename(new_dir_with_orig+".dbf", new_dir_with_code+".dbf")
            os.rename(new_dir_with_orig+".prj", new_dir_with_code+".prj")
            os.rename(new_dir_with_orig+".shp", new_dir_with_code+".shp")
            os.rename(new_dir_with_orig+".shx", new_dir_with_code+".shx")

            # Delete original zip file
            os.remove(zip_path)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def fetch_dem(shapefile=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset="30m", txt_file="download_urls.txt", save_to_txt=True, download_folder="dem_tiles", download=False, verbose=False):
    """
    Queries USGS API for DEM data (URLs) given specified parameters.

    This function targets USGS National Map API to fetch DEM data URLs using specified parameters and can either 
    save the list of URLs to a text file or download from the list of URLs immediately.

    Parameters
    ----------
    shapefile : str, optional
        Code of shapefile with which a bounding box will be generated (default is None). Overrides the 'bbox' parameter if set.
    bbox : dict, optional
        Bounding box coordinates to query for DEM data (default is {"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}).
    dataset : str, optional
        Code for DEM data to download (default is 30m).
    txt_file : str, optional
        Name of text file to save URLs to (default is download_urls.txt).
    save_to_txt : bool, optional
        Allow DEM URLs to be saved to a text file (default is True).
    download_folder : str, optional
        Name of the download folder (default is dem_tiles).
    download : bool, optional
        Allow DEM URLs retrieved to be immediately downloaded (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track download (default is False).

    Raises
    ------
    Exception
        If an invalid response while getting URLs off the USGS webpage, the function will error.
    """

    # Ensure dataset requested is valid
    if validate_codes([dataset], DATA_CODES) == -1: return
    
    # Get coordinate extents if a shape file was specified
    if shapefile is not None:
        if verbose is True: print("Reading in shapefile...")

        # Download shapefile
        download_shapefiles(shapefile)
        
        # Get path to shapefile
        shapefile_path = os.path.join(os.getcwd(), SHAPEFILE_FOLDER_NAME, shapefile, shapefile + ".shp")
        
        # Get extents of shape file
        coords = get_extents(shapefile_path)
        bbox["xmin"] = coords[0][0]
        bbox["ymax"] = coords[0][1]
        bbox["xmax"] = coords[1][0]
        bbox["ymin"] = coords[1][1]

    # Construct the query parameters
    if verbose is True: print("Setting boundary extents...")
    params = {
        "bbox": f"{bbox["xmin"]},{bbox["ymin"]},{bbox["xmax"]},{bbox["ymax"]}",
        "datasets": DATA_CODES[dataset],
        "prodFormats": DATA_FORMATS[dataset]
    }

    # Make a GET request to download DEM data
    if verbose is True: print("Requesting data from USGS...")
    base_url = "https://tnmaccess.nationalmap.gov/api/v1/products"
    response = requests.get(base_url, params=params)

    # Check for a successful request
    if response.status_code != 200:
        raise Exception(
            f"Failed to fetch data. Status code: {response.status_code}")

    # Convert JSON response to Python dictionary
    data = response.json()

    # Extract download URLs
    download_urls = [item["downloadURL"] for item in data["items"]]

    # Save URLs to text file 
    if save_to_txt is True:
        if verbose is True: print("Saving URLs to text file...")
        txt_path = os.path.join(os.getcwd(), txt_file) # Build full path to text file
        with open(txt_path, "w") as file:
            for url in download_urls:
                file.write(f"{url}\n")

    # Download the files from the URLs
    if download is True:
        download_files(download_urls, download_folder)

    # Completion message
    if verbose is True: print("Fetch process complete.")

####################################
### IMAGE MANIPULATION FUNCTIONS ###
####################################

def build_mosaic(input_folder, output_file, description=None, cleanup=False, verbose=False):
    """
    Builds a mosaic out of multiple GeoTIFF files.
    
    This function creates a mosaic from a list of GeoTIFF files utilizing the buildVRT() function from the GDAL library.

    Parameters
    ----------
    input_folder : str
        Name/path of the folder in data directory where files to mosaic together are located.
    output_file : str
        Name/path of mosaic file produced.
    description : str, optional
        Description to add to output raster band of mosaic file (default is None).
    cleanup : bool, optional
        Determines if files from `input_folder` should be deleted after mosaic is complete (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Create paths to files needed for computation
    vrt_path = os.path.join(os.getcwd(), VRT_DEFAULT_FILE_NAME)
    mosaic_path = determine_if_path(output_file)
    input_path = determine_if_path(input_folder)
    
    # Check for valid folder and put input files into list
    if verbose is True: print("Getting input files...")
    if validate_path_exists(input_path) == -1: return
    input_files = glob.glob(os.path.join(input_path, "*.tif"))

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
    os.remove(vrt_path)

    if verbose is True: print("Mosaic process complete.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def merge_shapefiles(input_folder, output_file, cleanup=False):
    """
    Merges shapefiles into a single shapefile.

    This function merges multiple shapefiles together into a single shapefile.
    Shapefiles provided should be .shp files.

    Parameters
    ----------
    input_folder : str
        Name of folder where shapefiles to merge are stored.
    output_file : str
        Name of output file that has merged shapefiles.
    cleanup : bool
        Determine if files from input folder should be deleted afterwards (default is False).
    """

    # Update path to input folder and output file if necessary
    input_path = determine_if_path(input_folder)
    output_path = determine_if_path(output_file)

    # Get all shapefiles
    input_files = glob.glob(os.path.join(input_path,"*.shp"))
    
    # Read and merge all shapefiles into a single GeoDataFrame
    merged_gdf = gpd.GeoDataFrame(pd.concat([gpd.read_file(file) for file in input_files], ignore_index=True))
    
    # Save the merged shapefile
    merged_gdf.to_file(output_file)

    # Cleanup input files
    if cleanup:
        shutil.rmtree(input_path)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def reproject(input_file, output_file, projection, cleanup=False, verbose=False):
    """
    Reprojects a GeoTIFF file to a specified projection.
    
    This function reprojects a specified GeoTIFF file to a new projection changing the coordinate system representation, 
    and the result is saved to a new file. Multithreading is utilized to improve performance.

    Parameters
    ----------
    input_file : str
        Name/path of GeoTIFF file in data directory to reproject.
    output_file : str
        Name/path of new reprojected GeoTIFF file to store in data directory.
    projection : str
        Projection to use for reprojection. Can be a EPSG code (e.g. EPSG:4326) or the path to a WKT file.
    cleanup : bool, optional
        Determines if `input_file` should be deleted after reprojection is complete (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Set full file paths for input and output
    input_path = determine_if_path(input_file)
    output_path = determine_if_path(output_file)

    # Ensure file to reproject exists
    if validate_path_exists(input_path) == -1: return
    
    if verbose is True: print("Reprojecting", os.path.basename(input_path) + "...")
    
    # Set warp options for reprojection and warp
    warp_options = gdal.WarpOptions(dstSRS=projection, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"],
                                    multithread=True, warpOptions=["NUM_THREADS=ALL_CPUS"])
    warp = gdal.Warp(output_path, input_path, options=warp_options)
    warp = None  # Close file

    # Delete files used for reprojection
    if cleanup is True:
        if verbose is True: print("Cleaning intermediary files...")
        os.remove(input_path)
        os.remove(input_path + ".aux.xml")

    # Completion message
    if verbose is True: print("Reprojection complete.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    input_path = determine_if_path(input_file)
    output_path = determine_if_path(output_file)
    
    # Window to crop by [left_x, top_y, width, height]
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

    # Perform translation
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def convert_file_format(input_path, output_path, new_format):
    """
    Converts a geospatial image into a specified format.

    Uses GDAL to convert a GDAL supported GeoTIFF file to SAGA supported SDAT, SGRD, and PRJ files or vice versa.

    Parameters
    ----------
    input_path : str
        Path to file to convert format to.
    output_path : str
        Path to store converted file in.
    new_format: str
        GDAL supported code with file type to convert file to.
    """

    # Check for valid format
    if new_format not in ["SAGA","GTiff"]:
        print("Invalid format passed. Use either 'SAGA' or 'GTiff'.")
        return
    
    # Assume format is SAGA by default, else update to GeoTIFF format options
    translate_options = gdal.TranslateOptions(format="SAGA")
    if new_format == "GTiff":
        translate_options = gdal.TranslateOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

    # Perform conversion
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_parameters_gdal(input_file, parameter_list):
    """
    Computes parameters using GDAL API.

    Computes parameters using GDAL API and updates the description
    of the final result to include the name of the parameter computed.

    Parameters
    ----------
    input_path : str
        Path to a GeoTIFF elevation file to compute parameters with.
    parameter_list : List[str]
        List of valid parameters to compute.
    """

    # Compute parameter for each one in list
    for param in parameter_list:
        output_file = os.path.join(os.path.dirname(input_file).replace('elevation_tiles',param+'_tiles'),os.path.basename(input_file))
        if param == 'aspect':
            dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            gdal.DEMProcessing(output_file, input_file, processing=param, options=dem_options)
        else:
            dem_options = gdal.DEMProcessingOptions(format='GTiff', creationOptions=['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'])
            gdal.DEMProcessing(output_file, input_file, processing=param, options=dem_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_and_compute_tile(window, items):
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
        the original raster file, the path to the cropped file, the buffer size, the
        method to compute the paramters with, and the list of parameters to compute,
        respectively.
    """
    
    # Sort out all items
    input_file = items[0]
    tile_file = items[1]
    buffer = items[2]
    method = items[3]
    convert = items[4]
    proj = items[5]
    params = items[6]

    # Crop tiles from input file and specified buffer and window
    window[0] = window[0] - buffer
    window[1] = window[1] - buffer
    window[2] = window[2] + (2*buffer)
    window[3] = window[3] + (2*buffer)
    crop_pixels(input_file, tile_file, window)

    # Compute parameters
    if method == 'SAGA':
        # Convert input tile to SGRD format then compute
        if convert:
            convert_file_format(tile_file, tile_file.replace('.tif','.sdat'), 'SAGA')
            tile_file = tile_file.replace('.tif','.sgrd')
        if 'all' in params:
            sw.compute_all_parameters(tile_file)
        else:
            sw.compute_parameters(tile_file, params)
    else:
        compute_parameters_gdal(tile_file, params)
        
    # Crop buffer region from computed tiles  
    params = SAGA_PARAMETER_LIST if 'all' in params else params
    for param in params:
        if (param != "channel_network") and (param != "drainage_basins"):
            # Set paths
            buffered_param_file = os.path.join(os.getcwd(),f"{param}_tiles",os.path.basename(tile_file))
            unbuffered_param_file = os.path.join(os.getcwd(),f"unbuffered_{param}_tiles",os.path.basename(tile_file.replace('.sgrd','.tif')))
            
            # Convert files back to GeoTIFF if needed
            if convert:
                convert_file_format(buffered_param_file.replace('.sgrd','.sdat'), buffered_param_file.replace('.sdat','.tif'), 'GTiff')
                buffered_param_file = buffered_param_file.replace('.sdat','.tif')
        
            # Crop off buffer
            ds = gdal.Open(buffered_param_file, 0)
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            ds = None
            param_window = [buffer, buffer, cols-(buffer*2), rows-(buffer*2)]
            crop_pixels(buffered_param_file, unbuffered_param_file, param_window)
        else:
            # Set path to shapefile
            param_file = os.path.join(os.getcwd(),f"{param}_tiles",os.path.basename(tile_file.replace('.tif','.shp').replace('.sgrd','.shp')))

            # Assign projection to shapefile
            gdf = gpd.read_file(param_file)
            gdf.set_crs(epsg=proj)
            gdf.to_file(param_file)
        
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_and_compute(input_file, column_length, row_length, parameter_list, compute_method='SAGA', buffer_size=10, num_processes=None, convert_file=True, projection=5070, cleanup=False, verbose=False):
    """
    Stages together important information to run cropping and computing.

    Computes desired windows for files to be cropped and computed as well as
    passes other important variables into another function that will do said
    functions. All windows and variables are passed into a Python multiprocessing
    pool, and will use a number of cores equivalent to the number of cropped files
    produced, or max cores if that number exceeds the number of cores.
    If the number of concurrent processes to use is not specified, the function
    will use the max number of cores available by default.

    Parameters
    ----------
    input_file : str
        Name or path to an input raster file to compute paramters from.
    column_length : int
        Column length of pixels to use for each cropped file.
    row_length : int
        Row length of pixels to use for each cropped file.
    parameter_list : List[str]
        List of paramters to compute. If special keyword 'all' in list, will compute parameters with SAGA using ta_compound 0.
    compute_method : str, optional
        API to use for computing terrain parameters (default is 'SAGA').
    buffer_size : int, optional
        Number of buffer pixels to use for cropping (default is 10).
    num_processes : int, optional
        Number of concurrent processes to use for computing (default is None).
    convert_file : bool, optional
        Determine if files should be converted to SGRD when using SAGA as the compute method (default is True).
    projection : int, optional
        Projection to set to metadata for shapefiles. Default is 5070.
    cleanup : bool, optional
        Specifies if cropped files used for computing parameters should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """

    # Get full path to input file if needed
    input_path = determine_if_path(input_file)
    
    # Create folders to store data intermediate data in
    input_tiles = os.path.join(os.getcwd(),'elevation_tiles')
    Path(input_tiles).mkdir(parents=True, exist_ok=True)
    if ('all' in parameter_list) and (compute_method == 'SAGA'):
        for parameter in SAGA_PARAMETER_LIST:
            Path(os.path.join(os.getcwd(),parameter+'_tiles')).mkdir(parents=True, exist_ok=True)
            if (parameter != "channel_network") and (parameter != "drainage_basins"):
                Path(os.path.join(os.getcwd(),'unbuffered_'+parameter+'_tiles')).mkdir(parents=True, exist_ok=True)
    else:
        for parameter in parameter_list:
            Path(os.path.join(os.getcwd(),parameter+'_tiles')).mkdir(parents=True, exist_ok=True)
            if (parameter != "channel_network") and (parameter != "drainage_basins"):
                Path(os.path.join(os.getcwd(),'unbuffered_'+parameter+'_tiles')).mkdir(parents=True, exist_ok=True)
    
    # Get number of columns and rows of data from input file
    ds = gdal.Open(input_path, 0)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    ds = None

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
        items.append((tile[0],[input_path,tile[1],buffer_size,compute_method,convert_file,projection,parameter_list]))

    # Setup multi-processing pool and compute
    if num_processes is None:
        num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(crop_and_compute_tile, items)

    # Cleanup cropped input tiles if specified
    if cleanup is True:
        if verbose is True: print("Cleaning intermediary files...")
        shutil.rmtree(input_tiles)

        parameter_list = SAGA_PARAMETER_LIST if ('all' in parameter_list) and (compute_method == 'SAGA') else parameter_list
        for parameter in parameter_list:
            if (parameter != "channel_network") and (parameter != "drainage_basins"):
                shutil.rmtree(os.path.join(os.getcwd(),f"{parameter}_tiles"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_region(codes, input_file, output_file):
    """
    Crops a raster file based off the region defined by a combination of shapefiles.
    
    This function uses a GDAL function to crop a raster file according to boundaries specified by a shapefile.
    Be careful that the bounds of the shapefile don't exceed the bounds of the data being cropped.
    The function automatically handles downloading shapefiles aren't present.

    Parameters
    ----------
    codes : List[str]
        List of shapefile codes that will outline the cropping region.
    input_file : str
        Name/path of input raster file to crop.
    output_file : str
        Name/path of cropped output raster file.
    """
    
    # Update file names to full paths
    input_path = determine_if_path(input_file)
    output_path = determine_if_path(output_file)

    # If codes are comma-separated, turn into a list
    if not isinstance(codes, list): 
        codes = list(codes.split(","))
    
    # Download shape files if any are missing
    download_shapefiles(codes)

    # Update path to shape files
    shape_paths = []
    for code in codes:
        shape_paths.append(os.path.join(os.getcwd(), SHAPEFILE_FOLDER_NAME, code, code + ".shp"))

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
    warp_options = gdal.WarpOptions(cutlineDSName=temp_combined_shp, cropToCutline=True, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
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

def crop_to_valid_data(input_file, output_file, projection='EPSG:4326', block_size=512):
    """
    Crops a border region of NaN values from a GeoTIFF file.

    This function uses blocking to scan through a GeoTIFF file to determine the extent of valid data 
    and crops excess Nan values at the borders of the spatial data.
    Blocking is used to help minimize RAM usage, and should be adjusted as accordingly.

    Parameters 
    ----------
    input_file : str
        Name/path of file to crop.
    output_file : str
        Name/path of file to save cropped data to.
    projection : str
        Name of projection to reference for translation of new file (default is EPSG:4326).
    block_size : int, optional
        Block size to use when computing extents (default is 512).
    """
    
    # Update file names to full paths
    input_path = determine_if_path(input_file)
    output_path = determine_if_path(output_file)
    
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
    
    out_ds = gdal.Translate(output_path, src_ds, projWin=[min_x, max_y, max_x, min_y], projWinSRS=projection)
    
    out_ds = None
    src_ds = None

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_coord(input_file, output_file, upper_left, lower_right):
    """
    Crops raster file to a specific region given specified coordinates.

    This function uses GDAL functions to crop a raster file based on a specified upper-left and lower-right coordinate.
    The 'upper_left' and 'lower_right' coordinates define the bounding box for cropping.
    Coordinates must be in the same projection as the raster.

    Parameters
    ----------
    input_file : str
        Name/path of the input raster file to crop.
    output_file : str
        Name of cropped raster to save.
    upper_left : tuple of float
       Float tuple specifying upper-left (x,y) coordinates to crop raster from.
    lower_right : tuple of float
        Float tuple specifying lower-right (x,y) coordinates to crop raster from.
    """
    
    # Build full paths if needed
    input_path = determine_if_path(input_file)
    output_path = determine_if_path(output_file)
    
    # Crop
    window = upper_left + lower_right
    translate_options = gdal.TranslateOptions(projWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
    gdal.Translate(output_path, input_path, options=translate_options)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_stack(input_folder, output_file):
    """
    Stacks multiple GeoTIFF files into a single GeoTIFF with multiple bands.

    This function takes multiple GeoTIFF files and combines them into one file represented by different bands. 
    Useful when multiple datasets need to be represented as one.
    The band order will be based off the 'input_files' list order.

    Parameters
    ----------
    input_folder : str
        Name/path of folder to GeoTIFF files to be stacked together.
    output_file : str
        Name/path of file to store stacked files to.
    """
    
    # Build full paths if needed
    input_path = determine_if_path(input_folder)
    output_path = determine_if_path(output_file)

    # Get input files
    input_files = glob.glob(os.path.join(input_path, "*.tif"))

    # Stack the files together
    vrt_options = gdal.BuildVRTOptions(separate=True)
    vrt = gdal.BuildVRT("stack.vrt", input_files, options=vrt_options)
    translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
    gdal.Translate(output_path, vrt, options=translate_options)
    vrt = None  # closes file

#####################################
### IMAGE VISUALIZATION FUNCTIONS ###
#####################################

def generate_img(tif, cmap="inferno", dpi=150, downsample=1, clean=False, title=None, nancolor="white",
                 ztype="Z", zunit=None, xyunit=None, vmin=None, vmax=None, reproject_gcs=False, shp_files=None, 
                 crop_shp=False, bordercolor="black", borderlinewidth=1.5, saveDir = None, verbose=False):
    """
    Plots a GeoTIFF image using matplotlib.

    This function visualizes a GeoTIFF file with optional parameters to ensure data visualization is made easily discernable and customizable.
    Using 'shp_files' without setting 'crop_shp' will allow you to plot the outline of the shapefile without cropping anything.
    Graph currently only tested for visualization with Jupyter Notebooks.
    Alternative colormaps can be found in the matplotlib documentation.

    Parameters
    ----------
    tif : str
        Name/path of the GeoTIFF file in the working directory to plot.
    cmap : str, optional
        Colormap used for visualization. Default is 'inferno'.
    dpi : int, optional
        Resolution in dots per inch for the figure. Default is 150.
    downsample : int, optional
        Factor to downsample the image by. Default is 10.
    clean : bool, optional
        If True, no extra data will be shown besides the plot image. Default is False.
    title : str, optional
        Title for the plot. Default will display the projection name.
    nancolor : str, optional
        Color to use for NaN values. Default is 'white'.
    ztype : str, optional
        Data that is represented by the z-axis. Default is 'Z'.
    zunit : str, optional
        Units for the data values (z-axis). Default is None and inferred from spatial reference.
    xyunit : str, optional
        Units for the x and y axes. Default is None and inferred from spatial reference.
    vmin : float, optional
        Value of lower bound for coloring on plot. Will be whatever the min value of the data is if none is specified.
    vmax : float, optional
        Value of upper bound for coloring on plot. Will be whatever the max value of the data is if none is specified.
    reproject_gcs : bool, optional
        Reproject a given raster from a projected coordinate system (PCS) into a geographic coordinate system (GCS) ESPG:4269.
    shp_files : str list, optional
        Comma-seperated list of strings with shape file codes to use for cropping. Default is None.
    crop_shp : bool, optional
        Flag to indicate if the shapefiles should be used for cropping. Default is False.
    bordercolor : str, optional
        Color for the shapefile boundary. Default is "black".
    borderlinewidth : float, optional
        Line width for the shapefile boundary. Default is 1.5.
    saveDir : str, optional
        Directory to save image to.
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).

    Outputs
    -------
    Image
        Displays image visualizing inputed GeoTIFF data with specified parameters.
    
    Returns
    -------
    raster_array: np.ndarray
        The raster array used for visualization. 
    """
    
    # Initial setup
    tif_dir_changed = False

    # Ensure file to plot exists
    tif_path = determine_if_path(tif)
    if validate_path_exists(tif_path) == -1: return
    
    # Reproject raster into geographic coordinate system if needed
    if reproject_gcs:
        reproject_path = os.path.join(os.getcwd(), "vis.tif")
        reproject(tif, "vis.tif", "EPSG:4269")
        if crop_shp is False:
            crop_path = os.path.join(os.getcwd(), "vis_trim_crop.tif")
            if verbose is True: print("Cropping NaN values...")
            crop_to_valid_data("vis.tif", "vis_trim_crop.tif")
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
            if verbose is True: print("Shapefile list is empty. Skipping shapefile cropping.")
        else:                
            if verbose is True: print("Cropping with combined shapefiles...")
            region_crop_path = os.path.join(os.getcwd(), "crop.tif")
            shape_paths = crop_region(shp_files, tif_path, "crop.tif")
            
            if tif_dir_changed:
                os.remove(tif_path)
            tif_path = region_crop_path
            tif_dir_changed = True

    if verbose is True: print("Reading in tif for visualization...")
    dataset = gdal.Open(tif_path)
    band = dataset.GetRasterBand(1)

    geotransform = dataset.GetGeoTransform()
    spatial_ref = osr.SpatialReference(wkt=dataset.GetProjection())

    # Extract spatial information about raster
    proj_name = spatial_ref.GetAttrValue("PROJECTION")
    proj_name = proj_name if proj_name else "GCS, No Projection"
    data_unit = zunit or spatial_ref.GetLinearUnitsName()
    coord_unit = xyunit or spatial_ref.GetAngularUnitsName()
    z_type = ztype if band.GetDescription() == "" else band.GetDescription()

    # Output verbose message if specified
    if verbose is True: 
        print(f"Geotransform:\n{geotransform}\n\nSpatial Reference:\n{spatial_ref}\n\nDocumentation on spatial reference format: https://docs.ogc.org/is/18-010r11/18-010r11.pdf\n")

    raster_array = gdal.Warp("", tif_path, format="MEM", 
                             width=int(dataset.RasterXSize/downsample), 
                             height=int(dataset.RasterYSize/downsample)).ReadAsArray()

    # Mask nodata values
    raster_array = np.ma.array(raster_array, mask=np.equal(raster_array, band.GetNoDataValue()))

    if verbose is True: print("Plotting data...")

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
        ax.axis("off")
    else:
        # Adjust colorbar and title
        cbar = fig.colorbar(sm, fraction=0.046*raster_array.shape[0]/raster_array.shape[1], pad=0.04)
        cbar_ticks = np.linspace(np.nanmin(raster_array), np.nanmax(raster_array), 8)
        cbar.set_ticks(cbar_ticks)
        cbar.set_label(f"{z_type} ({data_unit}s)")

        ax.set_title(title if title else f"Visualization of GEOTiff data using {proj_name}.", fontweight="bold")
        ax.tick_params(axis="both", which="both", bottom=True, top=False, left=True, right=False, color="black", length=5, width=1)

        ax.set_title(title or f"Visualization of GEOTiff data using {proj_name}.", fontweight="bold")

    # Set up the ticks for x and y axis
    x_ticks = np.linspace(ulx, lrx, 5)
    y_ticks = np.linspace(lry, uly, 5)

    # Format the tick labels to two decimal places
    x_tick_labels = [f"{tick:.2f}" for tick in x_ticks]
    y_tick_labels = [f"{tick:.2f}" for tick in y_ticks]

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_yticklabels(y_tick_labels)

    # Determine x and y labels based on whether data is lat-long or projected
    y_label = f"Latitude ({coord_unit}s)" if spatial_ref.EPSGTreatsAsLatLong() else f"Northing ({coord_unit}s)"
    x_label = f"Longitude ({coord_unit}s)" if spatial_ref.EPSGTreatsAsLatLong() else f"Easting ({coord_unit}s)"
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)

    # ax.ticklabel_format(style="plain", axis="both")  # Prevent scientific notation on tick labels
    ax.set_aspect("equal")

    if shp_files:
        for shp_file in shape_paths:
            overlay = gpd.read_file(shp_file)
            overlay.plot(ax=ax, edgecolor='black', facecolor='none')
            #overlay.boundary.plot(color=bordercolor, linewidth=borderlinewidth, ax=ax)

    if verbose is True: print("Done. Image should appear soon...")

    if saveDir is not None:
        fig.savefig(saveDir)

    if tif_dir_changed:
        os.remove(tif_path)

    return raster_array
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def visualize_shapefile(input_file, plot_title="", reproject_gcs=True, crop_to_shape=None):
    """
    Visualizes a shapefile.

    Visualizes a shapefile as a plot. It can be cropped to the bounds of another shapefile.

    Parameters
    ----------
    input_file : str
        Name of the shapefile in the working directory to plot.
    plot_title : str, optional
        Title to assign to plot (default is an empty string).
    reproject_crs : bool, optional
        Whether to reproject to standard EPSG:4269 coordinate system (default is True).
    crop_to_shape : str, optional
        State abbreviation to crop data to (default is None).

    Outputs
    -------
    Image
        Displays image visualizing inputed shapefile data with specified parameters.
    """
    
    # Load the input shapefile
    input_data = gpd.read_file(input_file)
    
    # Reproject shapefile
    if reproject_gcs:
        input_data = input_data.to_crs(epsg=4269) 

    # Crop data to another shapefile
    if crop_to_shape is not None:
        # Download shapefile if it does not exist
        download_shapefiles([crop_to_shape])
        
        crop_file = os.path.join(f"shapefiles/{crop_to_shape}/{crop_to_shape}.shp")
        crop_data = gpd.read_file(crop_file)

        # Ensure CRS match
        if input_data.crs != crop_data.crs:
            crop_data = crop_data.to_crs(epsg=4269)
        
        input_data = gpd.clip(input_data, crop_data)
    
    # Plot the shapefile(s)
    fig, ax = plt.subplots(
        figsize=(10, 10))
    input_data.plot(ax=ax)
    if crop_to_shape is not None:
        crop_data.plot(ax=ax, edgecolor='black', facecolor='none')
        
    plt.title(plot_title, fontsize=16, fontweight='bold')
    plt.xlabel("Longitude (Degrees)", fontsize=16)
    plt.ylabel("Latitude (Degrees)", fontsize=16)
    plt.xticks(fontsize=16) 
    plt.yticks(fontsize=16) 
    plt.show()

#############################
### FILE EXPORT FUNCTIONS ###
#############################

def extract_raster(csv_file, raster_file, band_names):
    """
    Extracts raster values and stores them with correlating coordinates found in a CSV file.

    This function reads raster values from a file and stores the value with correlating x and y coordinates in a CSV file.
    Order of 'band_names' should correlate to order of band names in specified raster file.
    If coordinates in CSV file are outside of bounds of raster, incorrect or no values will be extracted.

    Parameters
    ----------
    csv_file : str
        Name/path of CSV file to write to.
    raster_file : str
        Name/path of file to read raster values from.
    band_names : List[str]
        Names of bands to extract values from in input raster.
    """
    
    # Build full paths if needed
    csv_path = determine_if_path(csv_file)
    raster_path = determine_if_path(raster_file)
    
    # Extract values from raster corresponding to
    df = pd.read_csv(csv_path)

    ds = gdal.Open(raster_path, 0)
    gt = ds.GetGeoTransform()

    n_bands = ds.RasterCount
    bands = np.zeros((df.shape[0], n_bands))

    for i in range(df.shape[0]):
        px = int((df["x"][i] - gt[0]) / gt[1])
        py = int((df["y"][i] - gt[3]) / gt[5])

        for j in range(n_bands):
            band = ds.GetRasterBand(j + 1)
            val = band.ReadAsArray(px, py, 1, 1)
            bands[i, j] = val[0]
    ds = None

    for j in range(n_bands):
        df[band_names[j]] = bands[:, j]

    df.to_csv(csv_path, index=None)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def tif2csv(input_file, output_file="params", band_names=["elevation"]):
    """
    Converts a raster file into a CSV file.

    This function reads values from a GeoTIFF file and converts it into the CSV format. 
    The CSV columns are the x coordinate, y coordinate, and raster value.
    NaN values indicate no data found at a particular coordinate

    Parameters
    ----------
    input_file : str
        Name/path of GeoTIFF file to be stored as a CSV.
    output_file : str, optional
        Name/path of CSV file to save data to (default is 'params').
    band_names : List[str], optional
        Names of bands to pull data from (default is ['elevation']).
    """
    
    # Build full paths if needed
    input_path = determine_if_path(input_file)
    output_path = determine_if_path(output_file)

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

    column_names = ["x", "y"] + band_names
    stack = np.column_stack((x, y, bands))
    df = pd.DataFrame(stack, columns=column_names)
    df.dropna(inplace=True)
    df.to_csv(output_path, index=None)
