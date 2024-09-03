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
import tempfile
import zipfile
import shutil
import psutil
import math
import glob
import time
import json
import os
import re

# CONSTANTS
TEXT_FILE_EXTENSION = ".txt"
SAGA_FILE_EXTENSION = ".sdat"
GEOTIFF_FILE_EXTENSION = ".tif"
VRT_DEFAULT_FILE_NAME = "merged.vrt"
SHAPEFILE_FOLDER_NAME = "shapefiles"

# Silences a deprecation warning. 
gdal.UseExceptions()

###############################
### MISCELLANEOUS FUNCTIONS ###
###############################

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
    stdout, stderr = proc.communicate() # Get standard output and error of command

    # Print error message if execution returned error
    if proc.returncode != 0:
        raise RuntimeError("'%s' failed, error code: '%s', stdout: '%s', stderr: '%s'" % (
            ' '.join(arg_seq), proc.returncode, stdout.rstrip(), stderr.rstrip()))

###############################
### CONFIGURATION FUNCTIONS ###
###############################

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

def __validate_codes(input_codes, code_dictionary):
    for code in input_codes:
        if code not in list(code_dictionary.keys()):
            print("Code '{}' is not a valid code. Terminating execution.".format(code))
            return -1
    return 0

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __validate_path_exists(input_path):
    if not Path(input_path).exists():
        print("The path {} does not exist. Terminating execution.".format(input_path))
        return -1
    return 0

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __validate_process_count(core_count, file_count, process_count):
    if (process_count > core_count) or (process_count > file_count):
        print("The number of processes exceeds either the core count or input file count. Adjusting process count to lower of the two.")
        return min([core_count,file_count])
    return -1

#####################################
### FEARTURE EXTRACTION FUNCTIONS ###
#####################################

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

    # Get upper left and lower rights extents of shapefile and return
    ds = ogr.Open(shape_path)
    layer = ds.GetLayer()
    ext = layer.GetExtent()
    upper_left = (ext[0], ext[3])
    lower_right = (ext[1], ext[2])

    return upper_left, lower_right 

###############################
### FILE DOWNLOAD FUNCTIONS ###
###############################

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

    # Determine if URLs are a list or from a text file
    if type(download_list) is not list:
        download_list = os.path.join(os.getcwd(), download_list + TEXT_FILE_EXTENSION) # Update to be full path to text file
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
    Download specified shapefile(s) from a list of supplied codes.
    
    This function downloads specified shapefile(s) off the USGS webpage and stores them in a 'shapefiles' folder located in the data directory.
    All codes passed should be valid codes from the 'region_codes.txt' file.
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
    region_codes = __get_codes("region")
    if __validate_codes(codes, region_codes) == -1: return
    
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
            download_links.append(region_codes[code])
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

def fetch_dem(shapefile=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset="30m", txt_file="download_urls", save_to_txt=True, download_folder="dem_tiles", download=False, verbose=False):
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
        Name of text file to save URLs to (default is download_urls).
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
    data_codes = __get_codes("data")
    if __validate_codes([dataset], data_codes) == -1: return
    
    # Get coordinate extents if a shape file was specified
    if shapefile is not None:
        if verbose is True: print("Reading in shape file...")

        # Get extents of shape file
        coords = __get_extents(shapefile)
        bbox["xmin"] = coords[0][0]
        bbox["ymax"] = coords[0][1]
        bbox["xmax"] = coords[1][0]
        bbox["ymin"] = coords[1][1]

    # Construct the query parameters
    if verbose is True: print("Setting boundary extents...")
    params = {
        "bbox": f"{bbox["xmin"]},{bbox["ymin"]},{bbox["xmax"]},{bbox["ymax"]}",
        "datasets": data_codes[dataset],
        "prodFormats": "GeoTIFF"
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
        txt_path = os.path.join(os.getcwd(), txt_file + TEXT_FILE_EXTENSION) # Build full path to text file
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

def build_mosaic(input_folder, output_file, description, cleanup=False, verbose=False):
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
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Create paths to files needed for computation
    vrt_path = os.path.join(os.getcwd(), VRT_DEFAULT_FILE_NAME)
    mosaic_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    input_path = os.path.join(os.getcwd(), input_folder)
    
    # Check for valid folder and put input files into list
    if verbose is True: print("Getting input files...")
    if __validate_path_exists(input_path) == -1: return
    input_files = glob.glob(input_path + "/*" + GEOTIFF_FILE_EXTENSION)

    # Build VRT (mosaic files together)
    if verbose is True: print("Constructing VRT...")
    vrt = gdal.BuildVRT(vrt_path, input_files)
    translate_options = gdal.TranslateOptions(creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"])
    gdal.Translate(mosaic_path, vrt, options=translate_options)
    vrt = None  # close file

    # Update band description with name of terrain parameter
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

def reproject(input_file, output_file, projection, cleanup=False, verbose=False):
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
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Set full file paths for input and output
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)  

    # Ensure file to reproject exists
    if __validate_path_exists(input_path) == -1: return
    
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

def crop_into_tiles(input_file, output_folder, num_tiles, buffer=10, verbose=False):
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
    verbose : bool, optional
        Determine if additional print statements should be used to track computation (default is False).
    """
    
    # Update path to input file
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)

    # Ensure file to crop exists
    if __validate_path_exists(input_path) == -1: return
    
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

            tile_file = os.path.join(output_path, "tile_{0:04d}".format(tile_count) + GEOTIFF_FILE_EXTENSION)
            win = [j, i, ncols, nrows]

            # Upper left corner
            win[0] = max(0, win[0] - buffer)
            win[1] = max(0, win[1] - buffer)

            w = win[2] + 2 * buffer
            win[2] = w if win[0] + w < cols else cols - win[0]

            h = win[3] + 2 * buffer
            win[3] = h if win[1] + h < rows else rows - win[1]  

            __crop_pixels(input_file, os.path.join(output_folder, os.path.basename(tile_file).replace(GEOTIFF_FILE_EXTENSION, "")), win)
            if verbose is True: print(os.path.basename(tile_file), "cropped.")
            tile_count += 1 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __convert_file_format(input_path, output_path, new_format):
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
def __compute_params_saga(input_path, items):
    """
    Compute specified terrain parameters using SAGA API.

    This function computes specified terrain parameters using SAGA API within tmux sessions.
    Each parameter or parameters within a specified SAGA function are configured
    into its relevant bash command and run in a tmux session. Once the tmux 
    session is complete, it will move on to the next terrain parameter(s).
    The function handles conversion of the input elevation files from GeoTIFF to 
    SAGA supported SDAT and vice versa for computed parameters.
    
    Parameters
    ----------
    input_path : str
        Path to a GeoTIFF elevation file to compute parameters with.
    items : List[str]
        List of valid parameters to compute. The last four data in the list should be
        the number of cores to use for each SAGA instance, the output folder(s) prefix,
        and verbose option, respectively.
    """

    # Extract important variables from items list
    n_cores = items[-3]
    folder_prefix = items[-2]
    verbose = items[-1]
    param_list = items[:len(items)-3]

    # Elevation path config
    elev_path = os.path.dirname(input_path)
    elev_file_name = os.path.basename(input_path)
    elev_folder = elev_path[elev_path.rindex("/")+1:]
    saga_elev_path = elev_path.replace(elev_folder, "saga_" + elev_folder)
    Path(saga_elev_path).mkdir(parents=True, exist_ok=True)

    # Convert elevation file to SDAT
    if verbose is True: print("Converting", elev_file_name, "to SDAT...")
    saga_elev_file_name = elev_file_name.replace(GEOTIFF_FILE_EXTENSION, SAGA_FILE_EXTENSION)
    saga_input_path = os.path.join(saga_elev_path, saga_elev_file_name)
    if not Path(saga_input_path).exists():
        __convert_file_format(input_path, saga_input_path, "SAGA")
    
    # Create dictionaries of parameter directories and their associated parameter
    param_paths = {}
    saga_param_paths = {}
    for param in param_list:
        param_path = elev_path.replace(elev_folder, folder_prefix + "_" + param + "_tiles")
        Path(param_path).mkdir(parents=True, exist_ok=True)
        param_paths.update({param : param_path})
        saga_param_path = elev_path.replace(elev_folder, "saga_" + folder_prefix + "_" + param + "_tiles")
        Path(saga_param_path).mkdir(parents=True, exist_ok=True)
        saga_param_paths.update({param : saga_param_path})

    # Create dictionary of SAGA parameter paths with SGRD file extension for SAGA command line computation
    params_dict = {
        "elevation": saga_input_path.replace(SAGA_FILE_EXTENSION, ".sgrd")
    }

    # Add remaining parameters after elevation
    for param in param_list:
        params_dict.update({param: os.path.join(saga_param_paths[param], os.path.basename(saga_elev_file_name.replace(SAGA_FILE_EXTENSION, ".sgrd")))})

    
    # Start tracking execution time and memory usage
    file_name = os.path.basename(input_path)
    mem_log_file = os.path.join(os.getcwd(), 'mem_log.csv')
    f = open(mem_log_file, 'w')
    f.write('mem_usage\n')
    f.close()
    start_cmd = ['tmux', 'new-session', '-d', '-s', 'track-mem', 'python', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file]
    __bash(start_cmd)
    start_time = time.time()

    
    # Build base of command line function (tmux session and saga_cmd with core allocation call)
    cmd_base = ["saga_cmd", "-c="+str(n_cores)]
    
    # Slope, Aspect, and Curvature: https://saga-gis.sourceforge.io/saga_tool_doc/7.3.0/ta_morphometry_0.html
    if ("slope" in param_list) or ("aspect" in param_list) or ("profile_curvature" in param_list) or ("plan_curvature" in param_list):
        # Add command and elevation param
        cmd_curv = cmd_base + ["ta_morphometry", "0", "-ELEVATION", params_dict["elevation"]]

        # Add other requested parameters
        if "slope" in param_list:
            cmd_curv = cmd_curv + ["-SLOPE", params_dict["slope"]]
        if "aspect" in param_list:
            cmd_curv = cmd_curv + ["-ASPECT", params_dict["aspect"]]
        if "profile_curvature" in param_list:
            cmd_curv = cmd_curv + ["-C_PROF", params_dict["profile_curvature"]]
        if "plan_curvature" in param_list:
            cmd_curv = cmd_curv + ["-C_PLAN", params_dict["plan_curvature"]]

        if verbose is True: print("Computing curvature parameters...")
        __bash(cmd_curv) # Run
    
    if "hillshade" in param_list:
        # Build command
        cmd_shade = cmd_base + ["ta_lighting", "0", "-ELEVATION", params_dict["elevation"] , "-SHADE", params_dict["hillshade"]]

        if verbose is True: print("Computing analytical hillshading...")
        __bash(cmd_shade) # Run
    
    if "convergence_index" in param_list:
        # Build command
        cmd_ci = cmd_base + ["ta_morphometry", "1", "-ELEVATION", params_dict["elevation"] , "-RESULT", params_dict["convergence_index"]]

        if verbose is True: print("Computing convergence index...")
        __bash(cmd_ci) # Run

    # Convert SAGA files to GeoTIFF
    if verbose is True: print("Converting SAGA files to GeoTIFF...")
    for param in param_list:
        __convert_file_format(params_dict[param].replace(".sgrd",SAGA_FILE_EXTENSION), os.path.join(param_paths[param],elev_file_name), "GTiff")

    
    # End tracking
    ex_time = time.time() - start_time
    end_cmd = ['tmux','kill-session','-t','track-mem']
    __bash(end_cmd)

    # Compute peak memory usage and store results in a file
    mem_data = pd.read_csv(mem_log_file)
    df = pd.DataFrame(mem_data)
    mem_usage = df['mem_usage']
    peak_usage = max(mem_usage) - min(mem_usage)
    
    results_csv = os.path.join(os.getcwd(), 'individual_results.csv')
    results = [param_list[0], file_name, ex_time, peak_usage]
    formatted_result = param_list[0] + ',' + file_name + ',' + str(ex_time) + ',' + str(peak_usage) + '\n'
    f = open(results_csv, 'a')
    f.write(formatted_result)
    f.close()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_geotiled(input_folder, parameter_list, num_processes=4, num_cores=-1, output_folder_prefix='', cleanup=False, verbose=False):
    """
    Configures the multiprocessing pool for GEOtiled to begin computing terrain parameters.
    
    This function utilizes the multiprocessing library to allow for computation of parameters on different elevation GeoTIFF files at the same time.
    Parameters that can be computed are specified in the 'parameter_codes.txt', or the 'all' keyword can be passed to compute all parameters.
    It is better to keep `num_procs` as a low value for systems with low amounts of RAM.
    Note that multiprocessing isn't used when computing parameters with SAGA.
    
    Parameters
    ----------
    input_folder : str
        Name of the folder in the data directory containing DEM elevation files to compute terrain parameters from.
    param_list : List[str]
        List containing codes for terrain parameters to compute. The 'all' keyword will compute all terrain parameters.
    num_procs : int, optional
        Integer specifying the number of python instances to use for multiprocessing (default is 4).
    cores : int, optional
        Overrides automatic computation to determine number of cores to use for each SAGA computation (default is -1).
    output_folder_prefix : str, optional
        Prefix to attach to all output folders created for storing computed terrain paramters (default is '').
    cleanup : bool, optional
        Determine if elevation files should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """
    
    # Check to ensure input folder exists
    input_path = os.path.join(os.getcwd(), input_folder)
    if __validate_path_exists(input_path) == -1: return
    
    # Check for 'all' keyword in parameter list
    param_codes = __get_codes("parameter")
    if "all" in parameter_list:
        parameter_list = list(param_codes.keys())

    # Ensure all codes entered are valid codes and put full parameter names in different list
    param_names = []
    if __validate_codes(parameter_list, param_codes) == -1: return
    for param in parameter_list:
        param_names.append(param_codes[param])
    
    # Get files from input folder to compute parameters for 
    if verbose is True: print("Getting input files...")
    input_files = sorted(glob.glob(input_path + "/*" + GEOTIFF_FILE_EXTENSION))

    # Check to ensure number of processes doesn't exceed number of input files or cores
    new_num_procs = __validate_process_count(multiprocessing.cpu_count(), len(input_files), num_processes)
    if new_num_procs != -1: num_processes = new_num_procs

    # Determine how many cores each tile gets for computing parameters if not specified
    if num_cores == -1:
        num_cores = int(multiprocessing.cpu_count() / num_processes)
        if num_cores < 1: num_cores = 1

    # Combine all parameters and other variables into a single list to pass into multiprocessing pool
    item_info = param_names + [num_cores, output_folder_prefix, verbose]
    
    items = []
    for input_file in input_files:
        items.append((input_file, item_info))
    
    # Start the multiprocessing pool
    if verbose is True: print("Starting computation of parameters...")
    pool = multiprocessing.Pool(processes=num_processes) 
    pool.starmap(__compute_params_saga, items)

    # Remove files used to compute parameters
    if cleanup is True:
        if verbose is True: print("Cleaning files...")
        shutil.rmtree(input_path)
        # TODO - add all SDAT folders for removal

    # Successful completion message
    if verbose is True: print("GEOtiled computation done!")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def build_mosaic_buffer(input_folder, output_file, cleanup=False, verbose=False):
    """
    Builds mosaic from multiple GeoTIFF files that were cropped with a buffer region.

    This function is similar to the `build_mosaic` function but handles mosaicking together GeoTIFF files that were split to includes buffer regions
    by averaging points in the buffer regions together when merging.

    Parameters
    ----------
    input_folder : str
        Name of folder in data directory where files to mosaic together are located.
    output_file : str
        Name of mosaicked file produced.
    cleanup : bool, optional
        Determine if files used for mosaicking should be deleted after computation (default is False).
    verbose : bool, optional
        Determine if additional print statements should be used to track computation of parameters (default is False).
    """
    
    # Check to ensure input folder exists
    input_path = os.path.join(os.getcwd(), input_folder)
    if __validate_path_exists(input_path) == -1: return
    
    if verbose is True: print("Mosaicking process started...")
    
    # Get files from input folder to merge together
    input_files = glob.glob(input_path + "/*" + GEOTIFF_FILE_EXTENSION)

    # Build full path for VRT and output file
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)
    vrt_file = os.path.join(os.getcwd(), VRT_DEFAULT_FILE_NAME)

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
    __bash(cmd)

    # Remove intermediary files used to build mosaic
    if cleanup is True:
        if verbose is True: print("Cleaning files...")
        shutil.rmtree(input_path)
    os.remove(vrt_file)

    if verbose is True: print("Mosaic process completed.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

def __crop_region(codes, input_file, output_file):
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
        Name of input raster file to crop.
    output_file : str
        Name of cropped output raster file.
    """
    
    # Update file names to full paths
    input_path = os.path.join(os.getcwd(), input_file + GEOTIFF_FILE_EXTENSION)
    output_path = os.path.join(os.getcwd(), output_file + GEOTIFF_FILE_EXTENSION)

    # If codes are comma-separated, turn into a list
    if not isinstance(codes, list): 
        codes = list(codes.split(","))
    
    # Download shape files if any are missing
    __download_shapefiles(codes)

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

def __crop_to_valid_data(input_file, output_file, projection='EPSG:4326', block_size=512):
    """
    Crops a border region of NaN values from a GeoTIFF file.

    This function uses blocking to scan through a GeoTIFF file to determine the extent of valid data 
    and crops excess Nan values at the borders of the spatial data.
    Blocking is used to help minimize RAM usage, and should be adjusted as accordingly.

    Parameters 
    ----------
    input_file : str
        Name of file to crop.
    output_file : str
        Name of file to save cropped data to.
    projection : str
        Name of projection to reference for translation of new file (default is EPSG:4326).
    block_size : int, optional
        Block size to use when computing extents (default is 512).
    """
    
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
    
    out_ds = gdal.Translate(output_path, src_ds, projWin=[min_x, max_y, max_x, min_y], projWinSRS=projection)
    
    out_ds = None
    src_ds = None

#####################################
### IMAGE VISUALIZATION FUNCTIONS ###
#####################################

def generate_img(tif, cmap="inferno", dpi=150, downsample=1, clean=False, title=None, nancolor="green",
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
        Name of the GeoTIFF file in the working directory to plot.
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
        Color to use for NaN values. Default is 'green'.
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
        Reproject a given raster from a projected coordinate system (PCS) into a geographic coordinate system (GCS).
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
    tif_path = os.path.join(os.getcwd(), tif + GEOTIFF_FILE_EXTENSION)
    if __validate_path_exists(tif_path) == -1: return
    
    # Reproject raster into geographic coordinate system if needed
    if reproject_gcs:
        reproject_path = os.path.join(os.getcwd(), "vis" + GEOTIFF_FILE_EXTENSION)
        reproject(tif, "vis", "EPSG:4326")
        if crop_shp is False:
            crop_path = os.path.join(os.getcwd(), "vis_trim_crop" + GEOTIFF_FILE_EXTENSION)
            if verbose is True: print("Cropping NaN values...")
            __crop_to_valid_data("vis", "vis_trim_crop")
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
            region_crop_path = os.path.join(os.getcwd(), "crop" + GEOTIFF_FILE_EXTENSION)
            shape_paths = __crop_region(shp_files, os.path.basename(tif_path).replace(GEOTIFF_FILE_EXTENSION,""), "crop")
            
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
            overlay.boundary.plot(color=bordercolor, linewidth=borderlinewidth, ax=ax)

    if verbose is True: print("Done. Image should appear soon...")

    if saveDir is not None:
        fig.savefig(saveDir)

    if tif_dir_changed:
        os.remove(tif_path)

    return raster_array