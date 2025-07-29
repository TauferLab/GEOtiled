"""
GEOtiled Library v1.0.0
GCLab 2025

Compiled by Gabriel Laboy (@glaboy-vol) and Jay Ashworth (@washwor1)

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
import warnings
import requests
import zipfile
import shutil
import glob
import json
import math
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

# List of computable terrain parameters
COMPUTABLE_PARAMETERS = ["hillshade", "slope", "aspect", "plan_curvature", "profile_curvature", 
                         "convergence_index", "total_catchment_area", "flow_width", "specific_catchment_area", 
                         "channel_network", "channel_network_grid", "drainage_basins", "drainage_basins_grid", 
                         "flow_connectivity", "flow_direction", "channel_network_base_level", "channel_network_distance", 
                         "filled_depressions", "filled_flow_direction", "watershed_basins", "ls_factor", 
                         "topographic_wetness_index", "valley_depth", "relative_slope_position"]

# Downloaded datasets off the USGS
DATA_CODES = {"60m": "National Elevation Dataset (NED) Alaska 2 arc-second Current",
              "30m": "National Elevation Dataset (NED) 1 arc-second Current",
              "10m": "National Elevation Dataset (NED) 1/3 arc-second Current",
              "5m": "Alaska IFSAR 5 meter DEM",
              "3m": "National Elevation Dataset (NED) 1/9 arc-second",
              "1m": "Digital Elevation Model (DEM) 1 meter"}

# File formats of available datasets from the USGS
DATA_FORMATS = {"60m": "GeoTIFF",
                "30m": "GeoTIFF",
                "10m": "GeoTIFF",
                "5m": "ArcGrid",
                "3m": "IMG",
                "1m": "GeoTIFF"}

# Downloadable shapefile URLs from the USGS
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
### CONFIGURATION FUNCTIONS ###
###############################

def print_computable_parameters():
    """
    Outputs all computable parameters for GEOtiled from `COMPUTABLE_PARAMETERS`.
    A special indicator '(G)' indicates it is computable with GDAL.
    """

    try:
        print("Computable Terrain Parameters (G == GDAL only):")
        for parameter in COMPUTABLE_PARAMETERS:
            if (parameter in ["slope","aspect","hillshade"]):
                print(f"{parameter} (G)")
            else:
                print(parameter)
    except Exception as e:
        print(f"Error printing terrain parameters: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def extract_region_from_link(link):
    """
    Extracts the name of the region from USGS download link for shapefiles.

    Parameters
    ----------
    link : str
        USGS specific link that specifies URL to download a shapefile.

    Returns
    -------
    str
        Region name extracted from the link.
    """
    
    try:
        re_link = "https://prd-tnm.s3.amazonaws.com/StagedProducts/GovtUnit/Shape/GOVTUNIT_(.*)_State_Shape.zip"
        region = re.search(re_link, link).group(1).replace("_", " ") # Extract region from link
        return region
    except Exception as e:
        print(f"Error extracting region: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def print_region_codes():
    """
    Output codes and their correlating region from the `REGION_CODES` variable.
    These codes are meant to inform the user what shapefiles are available for download.
    """

    try:
        # Print all codes and their associated region
        print("Region codes and their associated region:")
        for key in REGION_CODES:
            region = extract_region_from_link(REGION_CODES[key])
            print(f"{key} : {region}")
    except Exception as e:
        print(f"Error printing codes: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def set_working_directory(path):
    """
    Creates needed directories then sets the specified path as the working directory.
    Important to run this before anything else to ensure data is stored in a desired location.

    Parameters
    ----------
    path : str
        Working directory to store data in.
    """

    try:
        # Create all needed directories first before setting the working directory
        os.makedirs(path, exist_ok=True)
        os.chdir(path) 
    except Exception as e:
        print(f"Error setting working directory: {e}")      
        return 

##################################
### INPUT VALIDATION FUNCTIONS ###
##################################

def validate_codes(codes, code_dictionary):
    """
    Validates if specified codes exist within correlating dictionary of values.
    
    Parameters
    ----------
    codes : List[str]
        List of codes to validate.
    code_dictionary : dict
        Contains full list of valid codes to compare against.
    """
    
    try:
        for code in codes:
            if code not in list(code_dictionary.keys()):
                print(f"{code} is not a valid code. Terminating execution.")
                exit()
    except Exception as e:
        print(f"Error validating codes: {e}")      
        return 

####################################
### FEATURE EXTRACTION FUNCTIONS ###
####################################

def get_shapefile_extents(shapefile):
    """
    Returns the bounding box (extents) of a shapefile.
    The extent is the upper left and lower right coordinates of the file.

    Parameters
    ----------
    shapefile : str
        Path to the shapefile to get extents for.
    
    Returns
    -------
    Tuples(float)
        Returns two tuples - the first being the upper left (x,y) coordinate, and the second being the lower right (x,y) coordinate
    """

    try:
        # Get upper left and lower rights extents of shapefile and return
        ds = ogr.Open(shapefile)
        layer = ds.GetLayer()
        ext = layer.GetExtent()
        upper_left = (ext[0], ext[3])
        lower_right = (ext[1], ext[2])
        ds = None # close file

        return upper_left, lower_right 
    except Exception as e:
        print(f"Error getting shapefile extents: {e}")      
        return 

###############################
### FILE DOWNLOAD FUNCTIONS ###
###############################

def get_file_size(url):
    """
    Retrieve the size of a file from a HEAD request to the provided URL and reads the 'Content-Length' header. 
    It's primarily designed to support the `download_files` function to calculate download sizes beforehand.

    Parameters
    ----------
    url : str
        A download URL.

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
    Download a single file from a URL.
    If the file already exists, skip the download.

    Parameters
    ----------
    url : str
        The download URL of the file.
    output_folder : str
        The path where the downloaded file will be stored.
    pbar : tqdm object
        Reference to the tqdm progress bar, used by the parent function to indicate download progress.

    Returns
    -------
    download_size : int
        The number of bytes downloaded.
    """
    
    try: 
        # If file already exists, download will be skipped, so total download size will be 0
        output_path = os.path.join(output_folder, url.split("/")[-1])
        if os.path.exists(output_path):
            return 0

        # Send request to download file and get download size
        response = requests.get(url, stream=True)
        downloaded_size = 0
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded_size += len(chunk)
                pbar.update(len(chunk))
                
        return downloaded_size
    except Exception as e:
        print(f"Error downloading file: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_files(download_list, download_folder):
    """
    Download file(s) from a list of URLs simultaneously using threading 
    and showcases download progress via a tqdm progress bar.

    Parameters
    ----------
    download_list : str, List[str]
        The name/path of a text file. This file should contain download URLs separated by newlines.
        A list of download URLs.
    download_folder : str
        Name/path of the folder where the downloaded files will be stored.
    """
    
    try:
        # Create the download folder
        Path(download_folder).mkdir(parents=True, exist_ok=True)
        
        # Determine if URLs are a list or from a text file
        if type(download_list) is not list:
            with open(download_list, "r", encoding="utf8") as dsvfile:
                urls = [url.strip().replace("'$", "") for url in dsvfile.readlines()]
        else:
            urls = download_list
        
        # Compute total size of files to download
        total_size = sum(get_file_size(url.strip()) for url in urls)

        # Create progress bar to track completion of downloads
        with tqdm(total=total_size, unit="B", unit_scale=True, ncols=100, desc="Downloading", colour="green") as pbar:
            with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                futures = [executor.submit(download_file, url, download_folder, pbar) for url in urls]
                for future in concurrent.futures.as_completed(futures):
                    size = future.result()
    except Exception as e:
        print(f"Error downloading files: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def download_shapefiles(codes):
    """
    Downloads specified shapefile(s) off the USGS webpage and stores them in a `shapefiles` folder located in the working directory.
    All codes passed should be valid codes from the REGION_CODES variable.
    It will skip downloading already existing shapefiles.

    Parameters
    ----------
    codes : str, List[str]
        Comma-separated string of valid codes.
        List of valid codes.
    """
    
    try:
        # If codes are a comma-separated string, turn it into a list
        if not isinstance(codes, list): 
            codes = list(codes.split(","))
        
        # Ensure codes passed are valid codes
        validate_codes(codes, REGION_CODES)

        # Iterate through all codes and download needed shapefiles
        download_links = []
        download_codes = []
        for code in codes:
            # Create path to shapefile
            shapefile_path = os.path.join(SHAPEFILE_FOLDER_NAME, code, f"{code}.shp")

            # Check if shapefile not already downloaded, and add download URL to list
            if not os.path.isfile(shapefile_path):
                # Get download URL based off region code
                download_links.append(REGION_CODES[code])
                download_codes.append(code)
    except Exception as e:
        print(f"Error fetching download URLs: {e}")      
        return 

    try:
        # Begin download
        if len(download_links) > 0:
            download_files(download_links, SHAPEFILE_FOLDER_NAME)

            # Unzip and rename files
            for iter, link in enumerate(download_links):
                # Isolate zip file name
                zip_name = os.path.basename(link)

                # Update path to zip file
                zip_path = os.path.join(SHAPEFILE_FOLDER_NAME, zip_name)

                # Unzip and get correct region files
                original_file_name = "GU_StateOrTerritory"
                original_file_path = f"Shape/{original_file_name}"
                with zipfile.ZipFile(zip_path, "r") as file:
                    file.extract(f"{original_file_path}.dbf", path=SHAPEFILE_FOLDER_NAME)
                    file.extract(f"{original_file_path}.prj", path=SHAPEFILE_FOLDER_NAME)
                    file.extract(f"{original_file_path}.shp", path=SHAPEFILE_FOLDER_NAME)
                    file.extract(f"{original_file_path}.shx", path=SHAPEFILE_FOLDER_NAME)
                
                file.close()
                
                # Rename folder and files
                dwnld_code = download_codes[iter]
                new_directory = os.path.join(SHAPEFILE_FOLDER_NAME, dwnld_code)
                new_dir_with_orig = os.path.join(new_directory, original_file_name)
                new_dir_with_code = os.path.join(new_directory, dwnld_code)
                os.rename(os.path.join(SHAPEFILE_FOLDER_NAME, "Shape"), new_directory)
                os.rename(f"{new_dir_with_orig}.dbf", f"{new_dir_with_code}.dbf")
                os.rename(f"{new_dir_with_orig}.prj", f"{new_dir_with_code}.prj")
                os.rename(f"{new_dir_with_orig}.shp", f"{new_dir_with_code}.shp")
                os.rename(f"{new_dir_with_orig}.shx", f"{new_dir_with_code}.shx")

                # Delete original zip file
                os.remove(zip_path)
    except Exception as e:
        print(f"Error downloading shapefiles: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def fetch_dems(shapefile=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset="30m", txt_file="download_urls.txt", save_to_txt=True, download_folder="dem_tiles", download=False, verbose=False):
    """
    Queries USGS National Map API to fetch DEM data URLs using specified filters and can either 
    save the list of URLs to a text file and/or download from the list of URLs immediately.

    Parameters
    ----------
    shapefile : str, optional
        Code of shapefile with which a bounding box will be generated (default is None). Overrides the `bbox` parameter if set.
    bbox : dict, optional
        Bounding box coordinates to query for DEM data (default is {"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}).
    dataset : str, optional
        Code for DEM data to download (default is '30m').
    txt_file : str, optional
        Name of text file to save URLs to (default is download_urls.txt).
    save_to_txt : bool, optional
        Allow DEM URLs to be saved to a text file (default is True).
    download_folder : str, optional
        Name of the download folder to store downloaded DEMs in (default is dem_tiles).
    download : bool, optional
        Allow DEM URLs retrieved to be immediately downloaded (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track download (default is False).
    """

    try:
        # Ensure dataset requested is valid
        if verbose: print("Starting fetch...")
        validate_codes([dataset], DATA_CODES)
        
        # Get coordinate extents if a shape file was specified
        if shapefile is not None:
            if verbose: print("Reading in shapefile...")

            # Download shapefile
            download_shapefiles(shapefile)
            
            # Get path to shapefile
            shapefile_path = os.path.join(SHAPEFILE_FOLDER_NAME, shapefile, f"{shapefile}.shp")
            
            # Get extents of shape file
            coords = get_shapefile_extents(shapefile_path)
            bbox["xmin"] = coords[0][0]
            bbox["ymax"] = coords[0][1]
            bbox["xmax"] = coords[1][0]
            bbox["ymin"] = coords[1][1]
    except Exception as e:
        print(f"Error reading shapefile: {e}")      
        return 

    try:
        # Construct the query parameters
        if verbose: print("Requesting data from USGS...")
        params = {
            "bbox": f"{bbox["xmin"]},{bbox["ymin"]},{bbox["xmax"]},{bbox["ymax"]}",
            "datasets": DATA_CODES[dataset],
            "prodFormats": DATA_FORMATS[dataset]
        }

        # Make a GET request to download DEM data
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
    except Exception as e:
        print(f"Error requesting data from the USGS: {e}")      
        return 

    try:
        # Save URLs to text file 
        if save_to_txt:
            if verbose: print("Saving URLs to a text file...")
            with open(txt_file, "w") as file:
                for url in download_urls:
                    file.write(f"{url}\n")
    except Exception as e:
        print(f"Error saving URLs to text file: {e}")      
        return 

    try:
        # Download the files from the URLs
        if download:
            if verbose: print("Downloading files...")
            download_files(download_urls, download_folder)
    except Exception as e:
        print(f"Error downloading files: {e}")      
        return 

####################################
### IMAGE MANIPULATION FUNCTIONS ###
####################################

def mosaic_rasters(input_folder, output_file, description=None, cleanup=False, verbose=False):
    """
    Mosaics together raster files into a single raster. The raster will only contain one band,
    and the user has the option to update the description of the band.

    Parameters
    ----------
    input_folder : str
        Name/path to the folder containing rasters to mosaic.
    output_file : str
        Name/path to the file to write the mosaicked rasters data to.
    description : str, optional
        Description to add to mosaicked raster's band (default is None).
    cleanup : bool, optional
        Determines if files from `input_folder` should be deleted after mosaic is complete (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    try:
        if verbose: print(f"Mosaicking files in {input_folder}...")
        
        # Check for valid folder and put input files into list
        input_files = glob.glob(os.path.join(input_folder, "*.tif"))

        # Create VRT to point to files to mosaic
        vrt = gdal.BuildVRT("merged.vrt", input_files)
        
        # Mosaic the files
        translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"])
        gdal.Translate(output_file, vrt, options=translate_options)
        vrt = None  # close file

        # Update band description
        if description is not None:
            if verbose: print("Updating band description...")
            dataset = gdal.Open(output_file)
            band = dataset.GetRasterBand(1)
            band.SetDescription(description)
            dataset = None  # close file

        # Delete intermediary tiles used to build mosaic
        if cleanup:
            if verbose: print("Cleaning intermediary files...")
            shutil.rmtree(input_folder)
        os.remove("merged.vrt")
    except Exception as e:
        print(f"Error mosaicking rasters: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def merge_shapefiles(input_folder, output_file, cleanup=False, verbose=False):
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
    cleanup : bool, optional
        Determine if files from input folder should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """

    try:
        if verbose: print("Merging shapefiles...")
        
        # Get all shapefiles
        input_files = glob.glob(os.path.join(input_folder,"*.shp"))
        
        # Read and merge all shapefiles into a single GeoDataFrame
        merged_gdf = gpd.GeoDataFrame(pd.concat([gpd.read_file(file) for file in input_files], ignore_index=True))
        
        # Save the merged shapefile
        merged_gdf.to_file(output_file)

        # Cleanup input files
        if cleanup:
            if verbose: print("Deleting input files...")
            shutil.rmtree(input_folder)
    except Exception as e:
        print(f"Error merging shapefiles: {e}")      
        return 
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def reproject(input_file, output_file, projection, cleanup=False, verbose=False):
    """
    Reprojects a raster file to a specified projection where the result is saved to a new file. 
    Multithreading is utilized to improve performance.

    Parameters
    ----------
    input_file : str
        Name/path of raster to reproject.
    output_file : str
        Name/path of file to write reprojected data to.
    projection : str
        Projection to use for reprojection. Can be a EPSG code (e.g. EPSG:4269) or the name/path to a WKT file.
    cleanup : bool, optional
        Determine if `input_file` should be deleted after reprojection is complete (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    try:
        if verbose: print(f"Reprojecting {input_file}...")
        
        # Set warp options for reprojection and warp
        warp_options = gdal.WarpOptions(dstSRS=projection, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"],
                                        multithread=True, warpOptions=["NUM_THREADS=ALL_CPUS"])
        warp = gdal.Warp(output_file, input_file, options=warp_options)
        warp = None  # Close file

        # Delete files used for reprojection
        if cleanup:
            if verbose: print("Cleaning intermediary files...")
            os.remove(input_file)
            if os.path.exists(f"{input_file}.aux.xml"):
                os.remove(f"{input_file}.aux.xml")
    except Exception as e:
        print(f"Error reprojecting file: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_to_window(input_file, output_file, window):
    """
    Crops a raster file to a specific window where the window references 2D matrix indices.

    Parameters
    ----------
    input_file : str
        Name/path to input file to be cropped.
    output_file : str
        Name/path to file to write the cropped data.
    window : List[int], Tuple(int)
        List or tuple of format [left_x, top_y, width, height] where left_x and top_y are pixel coordinates of the upper-left corner 
        of the cropping window, and width and height specify the dimensions of the cropping window in pixels.
    """
    
    try:
        # Window to crop by [left_x, top_y, width, height]
        translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

        # Perform translation
        gdal.Translate(output_file, input_file, options=translate_options)
    except Exception as e:
        print(f"Error cropping file: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def convert_file_format(input_file, output_file, new_format):
    """
    Converts a GDAL-supported raster to either GeoTIFF or SAGA-supported SDAT.
    Input file should be a GeoTIFF or SDAT file.

    Parameters
    ----------
    input_file : str
        Name/path to file to convert format of.
    output_file : str
        Name/path to write converted file to.
    new_format: str
        GDAL supported code with file type to convert file to.
    """

    try:
        # Check for valid format
        if new_format not in ["SAGA","GTiff"]:
            print("Invalid format passed. Use either 'SAGA' or 'GTiff'.")
            return
        
        # Assume format is SAGA by default, else update to GeoTIFF format options
        translate_options = gdal.TranslateOptions(format="SAGA")
        if new_format == "GTiff":
            translate_options = gdal.TranslateOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])

        # Perform conversion
        gdal.Translate(output_file, input_file, options=translate_options)
    except Exception as e:
        print(f"Error converting file format: {e}")      
        return 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_parameters_gdal(input_file, parameter_list):
    """
    Computes terrain parameters using GDAL API.

    Parameters
    ----------
    input_file : str
        Name/path to a raster file with elevation data to compute terrain parameters from.
    parameter_list : List[str]
        List of valid terrain parameters to compute.
    """

    try:
        # Compute each terrain parameter
        for parameter in parameter_list:
            output_file = os.path.join(os.path.dirname(input_file).replace("elevation_tiles", f"{parameter}_tiles"), os.path.basename(input_file))
            if parameter == "aspect":
                dem_options = gdal.DEMProcessingOptions(zeroForFlat=False, format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
                gdal.DEMProcessing(output_file, input_file, processing=parameter, options=dem_options)
            else:
                dem_options = gdal.DEMProcessingOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
                gdal.DEMProcessing(output_file, input_file, processing=parameter, options=dem_options)
    except Exception as e:
        print(f"Error computing terrain parameters with GDAL: {e}")      
        return  

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_and_compute_tile(window, items):
    """
    Designed to be used with the `crop_and_compute()` function.
    Crops a tile from a larger raster file with a specified window.
    Desired terrain parameters are then computed from the cropped tile.
    Debuffering of computed tiles is also handled.
    
    Parameters
    ----------
    window : List[int]
        Coordinates specifying what part of the original raster data to extract.
    items : List[]
        List of important variables passed from multiprocessing pool. Items include
        the input raster file, the name of the cropped tile file, the buffer size, the
        method to compute terrain paramters with, whether to convert files to SGRD, what
        projection to update shapefile metadata with, and the list of terrain parameters
        to compute, respectively.
    """
    
    try:
        # Sort out all items
        input_file = items[0]
        tile_file = items[1]
        buffer = items[2]
        method = items[3]
        convert = items[4]
        proj = items[5]
        params = items[6]
    except Exception as e:
        print(f"Error loading items: {e}")      
        return  

    try:
        # Crop tiles from input file and specified buffer and window
        window[0] = window[0] - buffer
        window[1] = window[1] - buffer
        window[2] = window[2] + (2*buffer)
        window[3] = window[3] + (2*buffer)
        crop_to_window(input_file, tile_file, window)
    except Exception as e:
        print(f"Error cropping tile: {e}")      
        return  

    try:
        # Compute terrain parameters using desired API
        if method == "SAGA":
            # Convert input tile to SGRD format
            if convert:
                convert_file_format(tile_file, tile_file.replace(".tif",".sdat"), "SAGA")
                tile_file = tile_file.replace(".tif",".sgrd")
            
            sw.compute_parameters(tile_file, params)
        else:
            compute_parameters_gdal(tile_file, params)
    except Exception as e:
        print(f"Error computing terrain parameters: {e}")      
        return  
        
    try:
        # Crop buffer region from computed tiles  
        for param in params:
            if (param != "channel_network") and (param != "drainage_basins"):
                # Set paths
                buffered_param_file = os.path.join(f"{param}_tiles",os.path.basename(tile_file.replace(".sgrd",".tif")))
                unbuffered_param_file = os.path.join(f"unbuffered_{param}_tiles",os.path.basename(tile_file.replace(".sgrd",".tif")))
                
                # Convert files back to GeoTIFF if needed
                if convert:
                    convert_file_format(buffered_param_file.replace(".tif",".sdat"), buffered_param_file, "GTiff")
            
                # Crop off buffer
                ds = gdal.Open(buffered_param_file, 0)
                cols = ds.RasterXSize
                rows = ds.RasterYSize
                ds = None
                param_window = [buffer, buffer, cols-(buffer*2), rows-(buffer*2)]
                crop_to_window(buffered_param_file, unbuffered_param_file, param_window)
            else:
                # Set path to shapefile
                param_file = os.path.join(f"{param}_tiles",os.path.basename(tile_file.replace(".tif",".shp").replace(".sgrd",".shp")))

                # Assign projection to shapefile
                gdf = gpd.read_file(param_file)
                gdf.set_crs(epsg=proj)
                gdf.to_file(param_file)
    except Exception as e:
        print(f"Error debuffering and assigning projections: {e}")      
        return  
        
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_and_compute(input_file, parameter_list, tile_dimensions=None, num_tiles=None, compute_method="SAGA", convert_file=True, projection=5070, buffer_size=10, num_processes=None, cleanup=False, verbose=False):
    """
    Computes terrain parameters in parallel. This function handles folder creation, cropping of the elevation data into tiles,
    parallel computation of terrain parameters, and debuffering of computed terrain parameter tiles.
    Concurrent computation will, by default, use the max number of available cores for processing 
    unless `num_processes` is specified.

    Parameters
    ----------
    input_file : str
        Name/path to input raster file with elevation data to compute terrain paramters from.
    parameter_list : List[str]
        List of terrain parameters to compute. If special keyword 'all' in list, will compute all terrain parameters for specific library.
    tile_dimensions : List[int], optional
        Column and row length of pixels to use for cropped tiles [x,y] (default is None).
    num_tiles : int, optional
        Will dynamically determine the size of the tiles if `column_length` and `row_length` is not specified (default is None).
    compute_method : str, optional
        API to use for computing terrain parameters (default is 'SAGA').
    convert_file : bool, optional
        Determine if files should be converted to SGRD when using SAGA as the compute method (default is True).
    projection : int, optional
        EPSG projection to set to metadata for created shapefiles (default is 5070).
    buffer_size : int, optional
        Number of buffer pixels to use for cropping (default is 10).
    num_processes : int, optional
        Number of concurrent processes to use for computing terrain parameters (default is None).
    cleanup : bool, optional
        Specifies if cropped files used for computing parameters should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """

    try:
        # Ensure user specified either a row and column length for tiles or a number of tiles
        if verbose: print("Creating folders...")
        if (tile_dimensions is None) and (num_tiles is None):
            print("Tile dimensions [column,row] or number of tiles must be specified")
            return

        # Create folders to store data intermediate data in
        Path("elevation_tiles").mkdir(parents=True, exist_ok=True)
        
        # Update parameter list if 'all' keyword specified
        if ("all" in parameter_list):
            if (compute_method == "SAGA"):
                parameter_list = COMPUTABLE_PARAMETERS
            elif (compute_method == "GDAL"):
                parameter_list = ["slope","aspect","hillshade"]
        
        # Create folders for all terrain parameters 
        for parameter in parameter_list:
            Path(os.path.join(f"{parameter}_tiles")).mkdir(parents=True, exist_ok=True)
            
            if (parameter != "channel_network") and (parameter != "drainage_basins"):
                Path(os.path.join(f"unbuffered_{parameter}_tiles")).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Error creating folders: {e}")      
        return  
    
    try:
        # Get number of columns and rows of data from input file
        if verbose: print("Computing windows...")
        ds = gdal.Open(input_file, gdal.GA_ReadOnly)
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        ds = None
        
        # Determine dimensions to use for tiles
        row_length, column_length = 1, 1
        if tile_dimensions is None:
            col_div = 1
            row_div = num_tiles
            diff = num_tiles - 1
            
            # Find product of tile count where difference between two values is minimized
            for divisor in range(2,int((num_tiles/2)+1),1):
                if num_tiles % divisor == 0:
                    remainder = num_tiles / divisor
                    temp_diff = abs(remainder - divisor)
                    if temp_diff < diff:
                        diff = temp_diff
                        col_div = divisor
                        row_div = remainder
                        
                        row_length = math.ceil(rows/row_div)            
                        column_length = math.ceil(cols/col_div)    
        else:
            row_length = tile_dimensions[1]
            column_length = tile_dimensions[0]

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
                tile_file = os.path.join("elevation_tiles", "tile_{0:04d}.tif".format(tile_count))
                tile_info.append([window, tile_file])
                tile_count += 1 
    except Exception as e:
        print(f"Error computing windows: {e}")   
        return

    try:
        # If compute method is GDAL, ensure file conversion is set to false
        if verbose: print("Starting computation...")
        if compute_method == "GDAL":
            convert_file = False

        # Store all required variables for multiprocessing into list
        items = []
        for tile in tile_info:
            items.append((tile[0],[input_file,tile[1],buffer_size,compute_method,convert_file,projection,parameter_list]))

        # Setup multi-processing pool and compute
        if num_processes is None:
            num_processes = len(items) if len(items) < os.cpu_count() else os.cpu_count()
        pool = multiprocessing.Pool(processes=num_processes)
        pool.starmap(crop_and_compute_tile, items)
    except Exception as e:
        print(f"Error computing terrain parameters: {e}")   
        return

    try:
        # Cleanup cropped input tiles if specified
        if cleanup:
            if verbose: print("Cleaning intermediary files...")
            shutil.rmtree("elevation_tiles")

            for parameter in parameter_list:
                if (parameter != "channel_network") and (parameter != "drainage_basins"):
                    shutil.rmtree(f"{parameter}_tiles")
                if (parameter == "watershed_basins"):
                    shutil.rmtree("pre_watershed_basins")
    except Exception as e:
        print(f"Error cleaning files: {e}")   
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_to_region(input_file, output_file, codes):
    """
    Crop a raster file to the boundaries specified by multiple shapefile.
    Multiple shapefiles can be given and will be combined into a single geometry for cropping.
    The function automatically handles downloading shapefiles aren't present.

    Parameters
    ----------
    input_file : str
        Name/path of input raster file to crop.
    output_file : str
        Name/path of output raster file cropped data will be written to.
    codes : List[str]
        List of shapefile codes that will outline the cropping region.
        
    Returns
    -------
    shape_paths : List[str]
        Returns a list of paths to shapefiles used for cropping.
    """
    
    try:
        # Download shape files if any are missing
        download_shapefiles(codes)
    except Exception as e:
        print(f"Error downloading shapefiles: {e}")   
        return

    try:
        # Update path to shape files
        shape_paths = []
        for code in codes:
            shape_paths.append(os.path.join(SHAPEFILE_FOLDER_NAME, code, f"{code}.shp"))

        # Read each shapefile, clean any invalid geometries, and union them
        gdfs = [gpd.read_file(shp_file).buffer(0) for shp_file in shape_paths]
        combined_geom = gdfs[0].unary_union
        for gdf in gdfs[1:]:
            combined_geom = combined_geom.union(gdf.unary_union)
        
        combined_gdf = gpd.GeoDataFrame(geometry=[combined_geom], crs=gdfs[0].crs)
            
        # Save the combined shapefile temporarily for cropping
        combined_gdf.to_file("temp_combined.shp")
    except Exception as e:
        print(f"Error combining geometries: {e}")   
        return

    try:
        # Do the cropping
        warp_options = gdal.WarpOptions(cutlineDSName="temp_combined.shp", cropToCutline=True, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
        warp = gdal.Warp(output_file, input_file, options=warp_options)
        warp = None
    except Exception as e:
        print(f"Error cropping file: {e}")   
        return

    try:
        # Remove the temporary combined shapefile
        os.remove("temp_combined.shp")
        os.remove("temp_combined.cpg")
        os.remove("temp_combined.dbf")
        os.remove("temp_combined.prj")
        os.remove("temp_combined.shx")

        # Return shapefile paths
        return shape_paths
    except Exception as e:
        print(f"Error removing temporary files: {e}")   
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_to_valid_data(input_file, output_file, projection='EPSG:4269', block_size=512):
    """
    Crops a raster file to the extents of valid data, where valid data is a row or column
    of data that contains at least one non-nan value (i.e., it crops borders of the raster file
    containing rows or columns of only nan values). Blocking is used to help minimize RAM usage, 
    and should be adjusted as accordingly.

    Parameters 
    ----------
    input_file : str
        Name/path of raster file to crop.
    output_file : str
        Name/path of raster file to save cropped data to.
    projection : str, optional
        Name of projection to reference for translation of new file (default is EPSG:4269).
    block_size : int, optional
        Block size to use when computing extents (default is 512).
    """
    
    try:
        # Open input file and get GeoTrasnform info
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
        
        # Crop file to coordinates 
        out_ds = gdal.Translate(output_file, src_ds, projWin=[min_x, max_y, max_x, min_y], projWinSRS=projection)
        
        # Close files
        out_ds = None
        src_ds = None
    except Exception as e:
        print(f"Error cropping raster to valid data: {e}")   
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def crop_to_coordinates(input_file, output_file, upper_left, lower_right):
    """
    Crops a raster file based on a specified upper-left and lower-right coordinates.
    The `upper_left` and `lower_right` coordinates define the bounding box for cropping.
    Coordinates must be in the same projection as the raster.

    Parameters
    ----------
    input_file : str
        Name/path of the input raster file to crop.
    output_file : str
        Name of cropped raster to save the cropped data to.
    upper_left : tuple of float
       Float tuple specifying upper-left (x,y) coordinates to crop raster from.
    lower_right : tuple of float
        Float tuple specifying lower-right (x,y) coordinates to crop raster from.
    """
    
    try:
        # Crop rastere to coordinates
        window = upper_left + lower_right
        translate_options = gdal.TranslateOptions(projWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
        gdal.Translate(output_file, input_file, options=translate_options)
    except Exception as e:
        print(f"Error cropping to coordinates: {e}")   
        return

#####################################
### IMAGE VISUALIZATION FUNCTIONS ###
#####################################

def plot_raster(input_file, plot_title=None, reproject_gcs=False, projection="EPSG:4269", shapefiles=None, remove_nans=False, 
                downsample=1, dpi=150, cmap="inferno", nancolor="white", ztype="Z", zunit=None, xyunit=None, vmin=None, 
                vmax=None, bordercolor="black", borderlinewidth=1.5, clean=False, save_to_img=None, verbose=False):
    """
    Plots a raster with user-specified modifications.

    Parameters
    ----------
    input_file : str
        Name/path of the raster file to plot.
    plot_title : str
        Title for plot (default is None).
    reproject_gcs : bool, optional
        Determine if raster should be reprojected before visualization (Default is False).
    projection : str, optional
        Projection to change raster to if reprojection is specified (default is 'EPSG:4269')
    shapefiles : List[str]
        List of valid shapefile codes to crop raster data to before visualization (default is None). 
    remove_nans : bool, optional
        Determine if nan-borders on data should be removed before plotting (default is False).
    downsample : int, optional
        Factor to downsample the raster by (default is 1).
    dpi : int, optional
        Resolution in dots per inch for the figure. Default is 150.
    cmap : str, optional
        Colormap used for visualization (default is 'inferno').
    nancolor : str, optional
        Color to use for NaN values (default is 'white').
    ztype : str, optional
        Data that is represented by the z-axis (default is 'Z').
    zunit : str, optional
        Units for the data values (z-axis) which is inferred unless specified (default is None).
    xyunit : str, optional
        Units for the x and y axes which is inferred unless specified (default is None).
    vmin : float, optional
        Value of lower bound for coloring on plot which is the minimum of the data unless specified (default is None).
    vmax : float, optional
        Value of upper bound for coloring on plot which is the maximum of the data unless specified (default is None).
    bordercolor : str, optional
        Color for the shapefile boundary if one is used for cropping (default is 'black').
    borderlinewidth : float, optional
        Line width for the shapefile boundary if one is used for cropping (default is 1.5).
    clean : bool, optional
        Determine whether to plot only the image (default is False).
    save_to_img : str, optional
        Name/path to image file to save plot to (default is None).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).

    Outputs
    -------
    Image
        Displays the plot of a raster with specified modifications.
    
    Returns
    -------
    raster_array : np.ndarray
        The data, in array format, of the plotted raster. 
    """
    
    try:
        # Reproject raster into geographic coordinate system 
        if reproject_gcs:
            if verbose: print("Reprojecting raster...")
            reproject(input_file, "temp_reprojected.tif", projection)
            input_file = "temp_reprojected.tif"
    except Exception as e:
        print(f"Error reprojecting raster: {e}")   
        return
    
    try:
        # Trim nan-border if desired
        if remove_nans:
            if verbose: print("Trimming NaN values...")
            crop_to_valid_data(input_file, "temp_trimmed.tif")
            input_file = "temp_trimmed.tif"
    except Exception as e:
        print(f"Error trimming nan values: {e}")   
        return

    try:
        # Crop to the desired shapefiles
        if shapefiles is not None:
            if verbose: print("Cropping to shapefiles...")
            shape_paths = crop_to_region(input_file, "temp_cropped.tif", shapefiles)
            input_file = "temp_cropped.tif"
    except Exception as e:
        print(f"Error cropping to shapefiles: {e}")   
        return

    try:
        # Begin reading in GeoTransform metadata for plotting
        if verbose: print("Reading raster metadata for visualization...")
        dataset = gdal.Open(input_file, gdal.GA_ReadOnly)
        band = dataset.GetRasterBand(1)

        geotransform = dataset.GetGeoTransform()
        spatial_ref = osr.SpatialReference(wkt=dataset.GetProjection())

        # Extract spatial information about raster
        proj_name = spatial_ref.GetAttrValue("PROJECTION")
        proj_name = proj_name if proj_name else "GCS, No Projection"
        data_unit = zunit or spatial_ref.GetLinearUnitsName()
        coord_unit = xyunit or spatial_ref.GetAngularUnitsName()
        z_type = ztype if band.GetDescription() == "" else band.GetDescription()

        # Print geostransform info if desired
        if verbose: 
            print(f"Geotransform:\n{geotransform}\n\nSpatial Reference:\n{spatial_ref}\n\nDocumentation on spatial reference format: https://docs.ogc.org/is/18-010r11/18-010r11.pdf\n")
    except Exception as e:
        print(f"Error loading GeoTransform information: {e}")   
        return

    try:
        # Downsample data
        if verbose: print("Downsampling data...")
        raster_array = gdal.Warp("", input_file, format="MEM", 
                                width=int(dataset.RasterXSize/downsample), 
                                height=int(dataset.RasterYSize/downsample)).ReadAsArray()

        # Mask nodata values
        raster_array = np.ma.array(raster_array, mask=np.equal(raster_array, band.GetNoDataValue()))
    except Exception as e:
        print(f"Error downsampling data: {e}")   
        return

    try:
        # Set up plotting environment
        if verbose: print("Starting plot construction...")
        cmap_instance = plt.get_cmap(cmap)
        cmap_instance.set_bad(color=nancolor)

        # Determine extent
        ulx, xres, _, uly, _, yres = geotransform
        lrx = ulx + (dataset.RasterXSize * xres)
        lry = uly + (dataset.RasterYSize * yres)

        # Plot
        fig, ax = plt.subplots(dpi=dpi)

        # Set min and max for plot colors
        vmn = vmin
        vmx = vmax
        if vmin is None:
            vmn = np.nanmin(raster_array)
        if vmax is None:
            vmx = np.nanmax(raster_array)
            
        # Show the plot
        sm = ax.imshow(raster_array, cmap=cmap_instance, vmin=vmn, vmax=vmx, extent=[ulx, lrx, lry, uly])
        
        # Determine if color bar and axis info should be shown
        if clean:
            ax.axis("off")
        else:
            # Adjust colorbar and title
            cbar = fig.colorbar(sm, fraction=0.046*raster_array.shape[0]/raster_array.shape[1], pad=0.04)
            cbar_ticks = np.linspace(np.nanmin(raster_array), np.nanmax(raster_array), 8)
            cbar.set_ticks(cbar_ticks)
            cbar.set_label(f"{z_type} ({data_unit}s)")

            ax.set_title(plot_title if plot_title else f"Visualization of GEOTiff data using {proj_name}.", fontweight="bold")
            ax.tick_params(axis="both", which="both", bottom=True, top=False, left=True, right=False, color="black", length=5, width=1)

        # Format tick labels
        x_ticks = np.linspace(ulx, lrx, 5)
        y_ticks = np.linspace(lry, uly, 5)
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
        ax.set_aspect("equal")

        # Plot shapefiles
        if shapefiles is not None:
            for shape_path in shape_paths:
                overlay = gpd.read_file(shape_path)
                overlay.plot(ax=ax, edgecolor=bordercolor, facecolor='none', linewidth=borderlinewidth)

        # Save plot to a file
        if save_to_img is not None:
            fig.savefig(save_to_img)
    except Exception as e:
        print(f"Error plotting data: {e}")   
        return

    try:
        # Cleanup temporary files
        if verbose: print("Cleaning intermediary files...")
        if reproject_gcs:
            os.remove("temp_reprojected.tif")
        if remove_nans:
            os.remove("temp_trimmed.tif")
        if shapefiles is not None:
            os.remove("temp_cropped.tif")

        # Return array data
        return raster_array
    except Exception as e:
        print(f"Error cleaning intermediary files: {e}")   
        return
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_shapefile(input_file, plot_title="", reproject_gcs=False, projection=4269, crop_to_shape=None, save_to_img=None, verbose=False):
    """
    Visualizes a shapefile in a plot. The projection can be changed, and the data can be cropped to the bounds of another shapefile.
    The user also has the option to save the plot to a file.

    Parameters
    ----------
    input_file : str
        Name/path to the shapefile to plot.
    plot_title : str, optional
        Title to assign to plot (default is '').
    reproject_gcs : bool, optional
        Determine whether to reproject data before visualization (default is False).
    projection : int, optional
        If reprojecting data, specifies which EPSG to reproject to (default is 4269).
    crop_to_shape : str, optional
        State shapefile to crop data to (default is None).
    save_to_img : str, optional
        Name/path of file to save plot to (default is None).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).

    Outputs
    -------
    Image
        Creates a plot for the shapefile.
    """
    
    try:
        # Load the shapefile
        if verbose: print("Reading in shapefile...")
        input_data = gpd.read_file(input_file)

        # Reproject shapefile
        if reproject_gcs:
            if verbose: print("Reprojecting...")
            input_data = input_data.to_crs(epsg=projection) 

        # Crop data to another shapefile
        if crop_to_shape is not None:
            # Download shapefile if not already done so
            if verbose: print("Cropping to another shapefile...")
            download_shapefiles([crop_to_shape])
            
            crop_file = os.path.join(f"shapefiles/{crop_to_shape}/{crop_to_shape}.shp")
            crop_data = gpd.read_file(crop_file)

            # Ensure CRS match
            if input_data.crs != crop_data.crs:
                crop_data = crop_data.to_crs(epsg=projection)
            
            input_data = gpd.clip(input_data, crop_data)
    except Exception as e:
        print(f"Error loading and modifying shapefile: {e}")
        return
    
    try:
        # Get bounds of shapefile
        if verbose: print("Constructing plot...")
        minx, miny, maxx, maxy = input_data.total_bounds
        
        # Plot the shapefile(s)
        fig, ax = plt.subplots(figsize=(10, 10))
        
        input_data.plot(ax=ax)
        if crop_to_shape is not None:
            crop_data.plot(ax=ax, edgecolor='black', facecolor='none')
            
            # Update bounds to cropped shape
            minx, miny, maxx, maxy = crop_data.total_bounds
            
        # Set up the ticks for x and y axis
        x_ticks = np.linspace(minx, maxx, 5)
        y_ticks = np.linspace(miny, maxy, 5)

        # Format the tick labels to two decimal places
        x_tick_labels = [f"{tick:.2f}" for tick in x_ticks]
        y_tick_labels = [f"{tick:.2f}" for tick in y_ticks]
            
        plt.title(plot_title, fontsize=16, fontweight='bold')
        plt.xlabel("Longitude (Degrees)", fontsize=16)
        plt.ylabel("Latitude (Degrees)", fontsize=16)
        ax.set_xlim(minx, maxx)
        ax.set_ylim(miny, maxy)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_xticklabels(x_tick_labels, fontsize=16)
        ax.set_yticklabels(y_tick_labels, fontsize=16)
        plt.show()
        
        if save_to_img is not None:
            plt.savefig(save_to_img)
    except Exception as e:
        print(f"Error plotting shapefile: {e}")
        return
    

#############################
### FILE EXPORT FUNCTIONS ###
#############################

def build_stack(input_files, output_file, verbose=False):
    """
    Stacks multiple rasters into a single raster with multiple bands.
    The band order will be based off the `input_files` list order.

    Parameters
    ----------
    input_files : List[str]
        Name/paths of rasters to stack together.
    output_file : str
        Name/path of file to store stacked raster data to.
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """
    
    try:
        # Create VRT file
        if verbose: print("Stacking rasters...")
        vrt_options = gdal.BuildVRTOptions(separate=True)
        vrt = gdal.BuildVRT("stack.vrt", input_files, options=vrt_options)
        
        # Stack files together
        translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
        gdal.Translate(output_file, vrt, options=translate_options)
        vrt = None  # closes file
        
        os.remove("stack.vrt")
    except Exception as e:
        print(f"Error stacking rasters: {e}")   
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def extract_raster(csv_file, raster_file, band_names, verbose=False):
    """
    Uploads raster data to a CSV that already has x,y coordinates. Only data from the raster that has matching 
    x,y coordinates are written to the CSV. Order of `band_names` should correlate to order of the raster's bands.
    If no x,y coordinates in the CSV correlate to those of the raster file, incorrect or no values will be extracted.

    Parameters
    ----------
    csv_file : str
        Name/path to CSV file to read/write to.
    raster_file : str
        Name/path to raster file to read to values from.
    band_names : List[str], optional
        Names of columns correlating to raster bands. Order matters.
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """
    
    try:
        # Load in CSV data to a dataframe
        if verbose: print("Loading in CSV...")
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error loading in CSV to dataframe: {e}")
        return

    try:
        # Load in raster data and get Geotransform info
        if verbose: print("Loading in raster data...")
        ds = gdal.Open(raster_file, 0)
        gt = ds.GetGeoTransform()

        n_bands = ds.RasterCount
        bands = np.zeros((df.shape[0], n_bands))

        # Get data from raster correlating to existing x,y coordinate from the CSV file
        if verbose: print("Formatting raster data...")
        for i in range(df.shape[0]):
            px = int((df["x"][i] - gt[0]) / gt[1])
            py = int((df["y"][i] - gt[3]) / gt[5])

            for j in range(n_bands):
                band = ds.GetRasterBand(j + 1)
                val = band.ReadAsArray(px, py, 1, 1)
                bands[i, j] = val[0]
        ds = None # close file
    except Exception as e:
        print(f"Error loading and formatting raster data: {e}")
        return

    try:
        # Format dataframe with new data and write back to CSV
        if verbose: print("Writing data to CSV...")
        for j in range(n_bands):
            df[band_names[j]] = bands[:, j]

        df.to_csv(csv_file, index=None)
    except Exception as e:
        print(f"Error writing raster results to CSV: {e}")
        return

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def tif2csv(input_file, output_file, band_names, verbose=False):
    """
    Converts a GeoTIFF file and converts it into the CSV format. The CSV columns are the x coordinate, y coordinate, and raster values.
    NaN values indicate no data found at a particular index. If the GeoTIFF has multiple bands, the user can set the names of subsequent 
    columns in the CSV containing each bands value for each x,y coordinate using the `band_names` parameter.

    Parameters
    ----------
    input_file : str
        Name/path to GeoTIFF file to convert to a CSV.
    output_file : str
        Name/path to CSV file to save CSV data to.
    band_names : List[str], optional
        Names of columns correlating to raster bands. Order matters.
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """
    
    try:
        # Open raster and extract geospatial coordinate data
        if verbose: print("Loading in raster data...")
        ds = gdal.Open(input_file, 0)
        xmin, res, _, ymax, _, _ = ds.GetGeoTransform()
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        xstart = xmin + res / 2
        ystart = ymax - res / 2

        # Convert x,y matrix indices to longitude,latitude coordinates
        x = np.arange(xstart, xstart + xsize * res, res)
        y = np.arange(ystart, ystart - ysize * res, -res)
        x = np.tile(x[:xsize], ysize)
        y = np.repeat(y[:ysize], xsize)

        # Get data values for each band
        if verbose: print("Formatting raster data...")
        n_bands = ds.RasterCount
        bands = np.zeros((x.shape[0], n_bands))
        for k in range(1, n_bands + 1):
            band = ds.GetRasterBand(k)
            data = band.ReadAsArray()
            data = np.ma.array(data, mask=np.equal(data, band.GetNoDataValue()))
            data = data.filled(np.nan)
            bands[:, k-1] = data.flatten()
        ds = None # close file
    except Exception as e:
        print(f"Error loading and formatting GeoTIFF data: {e}")
        return

    try:
        if verbose: print("Formatting and writing data to CSV...")
        # Format and load all values into a dataframe
        column_names = ["x", "y"] + band_names
        stack = np.column_stack((x, y, bands))
        df = pd.DataFrame(stack, columns=column_names)
        df.dropna(inplace=True)
        
        # Save dataframe to a CSV
        df.to_csv(output_file, index=None)
    except Exception as e:
        print(f"Error saving dataframe to CSV: {e}")
        return