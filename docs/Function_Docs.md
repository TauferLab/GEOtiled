# GEOtiled Function Documentation


## Dictionary Codes


### USGS_DATASET_CODES

Codes specifying which dataset to download GeoTIFF elevation files for off the [USGS webpage](https://apps.nationalmap.gov/downloader/#/).

Currently supported codes and their correlating dataset:
* 30m (National Elevation Dataset (NED) 1 arc-second Current)
* 10m (National Elevation Dataset (NED) 1/3 arc-second Current)


### TERRAIN_PARAMETER_CODES

Codes specifying terrain parameters to be computed by GEOtiled.

Currently supported codes and their correlating parameter:
* SLP (slope)
* ASP (aspect)
* HLSD (hillshade)


### SHAPE_FILE_CODES

Codes specifying a US state or territory for which to fetch the shape file for off the [USGS webpage](https://apps.nationalmap.gov/downloader/#/).

Codes specifying the region and their correlating parameter:
* AL (Alabama)
* AK (Alaska)
* AR (Arkansas)
* AZ (Arizona)
* CA (California)
* CO (Colorado)
* CT (Connecticut)
* DC (District of Columbia)
* DE (Delaware)
* FL (Florida)
* GA (Georgia)
* GU (Guam)
* HI (Hawaii)
* ID (Idahoo)
* IL (Illinois)
* IN (Indiana)
* IA (Iowa)
* KS (Kansas)
* KY (Kentucky)
* LA (Louisiana)
* MA (Massachusetts)
* MD (Maryland)
* ME (Maine)
* MI (Michigan)
* MN (Minnesota)
* MO (Missouri)
* MP (Northern Mariana Islands)
* MS (Mississippi)
* MT (Montana)
* NC (North Carolina)
* ND (North Dakota)
* NE (Nebraska)
* NH (New Hamshire)
* NJ (New Jersey)
* NM (New Mexico)
* NV (Nevada)
* NY (New York)
* OH (Ohio)
* OK (Oklahoma)
* OR (Oregon)
* PA (Pennsylvania)
* PR (Puerto Rico)
* RI (Rhode Island)
* SC (South Carolina)
* SD (South Dakota)
* TN (Tennessee)
* TX (Texas)
* UT (Utah)
* VA (Virginia)
* VI (Virgin Islands)
* VT (Vermont)
* WA (Washington)
* WI (Wisconsin)
* WV (West Virginia)
* WY (Wyoming)


## Functions


### `bash(argv)`

#### Executes a command in bash.
This function acts as a wrapper to execute bash commands using the subprocess.Popen() method. Commands are executed synchronously, and errors are caught and raised.

#### Required Parameters
argv : str list
* List of arguments for a bash command.
* The list should be ordered the same way you would write the function in a command line (e.g., ["ls", "-lh", "~/"]).

#### Outputs
Command Outputs
* Any output(s) produced by the bash command.

#### Error States
RuntimeError
* Will raise a RuntimeError if Popen() returns an error and print out the error code, stdout, and stderr.


### `build_mosaic(input_folder, output_file, description, cleanup=False)`

#### Builds a mosaic out of multiple GeoTIFF files.
This function creates a mosaic from a list of GeoTIFF files utilizing the `buildVRT()` function from the GDAL library.

#### Required Parameters
input_folder : str
* String specifying name of folder in data directory where GeoTIFF files to mosaic together are located.

output_file : str
* String specifying name of mosaicked file produced.

description : str
* String specifying description to add to output raster band of final GeoTIFF file.

#### Optional Parameters
cleanup : bool
* Boolean specifying if tiles in input_folder used for computation should be deleted after computation is complete. Default is False.

#### Outputs
GeoTIFF File
* Produces a GeoTIFF file that is the mosaic of all GeoTIFF files from the specified input folder.


### `compute_geotiled(input_folder, param_list, num_procs, cleanup=False)`

#### Configures the multiprocessing pool for GEOtiled to begin computing terrain parameters.
This function utilizes the multiprocessing library to allow for computation of parameters on different elevation GeoTIFF files at the same time.

#### Required Parameters


#### Optional Parameters


#### Notes
* It is best practice set 'num_procs' to a small value if the system does not have a lot of RAM.


### `compute_params(input_file, param_list)`

#### Computes terrain parameters for a given elevation GeoTIFF.
This function utilizes GDAL and GRASS libraries to compute parameters like slope, aspect, etc. from a provided parameter list.

#### Required Parameters
input_file : str
* Path to a GeoTIFF elevation file to compute params with.

param_list : str list
* String ist of valid parameters to compute for.

#### Outputs
Folder(s)
* Will create folders for computed params if they do not already exist. Folders following the naming scheme of 'parameter_tiles'.

File(s)
* The computed GeoTIFF files of all parameters specified stored in their appropiate folder.

#### Notes
* GDAL is used to compute slope, aspect, and hillshade.
* GRASS is used to compute all other parameters.


### `crops_into_tiles(input_file, output_folder, num_tiles, buffer=10)`

#### Splits a GeoTIFF file into smaller, equally-sized tiles.
This function divides a GeoTIFF file into a specified number of tiles with added buffer regions to assist with rapid computation of parameters by GEOtiled.

#### Required Parameters
input_file : str
* Name of the GeoTIFF file in the data directory to crop.

output_folder : str
* Name of the folder in the data directory to store the cropped tiles.

n_tiles : int
* Number of total tiles to produce. Should be a perfect square number.

#### Optional Parameters
* Specifies the buffer size (number of added pixels at the edges of the tiles generated. Default is 10.

#### Outputs
Folder
* Folder to hold cropped tiles in data directory where name is specified by 'output_folder' variable.
Cropped Files
* The cropped GeoTIFF files stored in the 'output_folder'.


### `crop_pixels(input_file, output_file, window)`

#### Crops a raster file to a specific region given a specified window.
This function uses GDAL functions to crop data based off pixel coordinates rather than geospatial coordinates.

#### Required Parameters
input_file : str
* String specifying path to input file to be cropped.

output_file : str
* String specifying path to new output file produced after cropping.

window : int list | int tuple
* List or tuple in the format [left_x, top_y, width, height] where left_x and top_y are pixel coordinates of the upper-left corner of the cropping window, and width and height specify the dimensions of the cropping window in pixels.

#### Outputs
Cropped File
* Cropped raster file saved at the specified output_file path.

#### Notes
* Ensure the specified pixel window is within bounds of the input raster or an error will be raised.


### `download_file(url, folder, pbar)`

#### Downloads a file found at a URL and stores it in the specified folder.
This function facilitates download of a single file and is used by the `download_files` function.

#### Required Parameters
url : str
* String specifying the URL to download the file from.

folder : str
* String specifying folder in data directory where file will be stored.

pbar : tqdm object
* Reference to a tqdm progress bar to indiciate download progress.

#### Outputs
File
* Downloaded file in the specified folder.

#### Returns
int
* Integer specifying the number of bytes downloaded.

#### Notes
* This function is designed to be used specifically by `download_files`
* If the file being downloaded already exists, no download occurs, and the function returns 0.


### `download_files(download_list, download_folder)`

#### Downloads file(s) from a list of provided URLs and stores them in a specified folder.
This function allows for simultaneous download of files off the [USGS webpage](https://apps.nationalmap.gov/downloader/#/) via threading and visualizes progress via a tqdm bar.

#### Required Parameters
download_list : str | str list
* A string specifying the name of a text file with newline separated download URLs.
* A list of string with download URLs.

download_folder : str
* String denoting folder in data directory where downloaded files will be stored.

#### Outputs
Folder
* Will create the specified `download_folder` if it doesn't already exist

#### Notes
* The function uses `ThreadPoolExecutor` from the `concurrent.futures` library to do multi-threaded downloads for efficiency.
* Will not download files if the file already exists, but the progress bar will not reflect it. 


### `download_shape_files(codes)`

#### Download shape file(s) from a list of specified codes.
This function downloads the specified shape file(s) off the [USGS webpage](https://apps.nationalmap.gov/downloader/#/) and stores them in a 'shape_files' folder located in the data directory.

#### Required Parameters
codes : str list
* String list with codes specified by the `SHAPE_FILE_CODES` variable.

#### Outputs
File(s)
* Downloaded shape files located in 'shape_files' folder.

#### Notes
* If the shape file already exists, it will not redownload it.
* The function uses the `download_files` function for dowloading.

#### Error States
* If the shape file is unable to be retrieved from the USGS webpage, a response error will be raised.


### `fetch_dem(shape_file=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset='30m', prod_format='GeoTIFF', txt_file='download_urls', save_to_txt=True, download=False)`

#### Queries USGS API for DEM data (URLs) given specified parameters.
This function targets USGS National Map API to fetch DEM data URLs using specified parameters and can either save the list of URLs to a text file or download from the list of URLs immediately.

#### Optional Parameters
shape_file : str list
* String list with single shapefile code to generate bounding box. Default is None.

bbox : float dictionary
* Dictionary containing bounding box coordinates to query for DEM data. Default is {"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}.

dataset : str
* String containing code for DEM data to download. Default is 30m.

prod_format : str
* File type of DEM data to download. Default is GeoTIFF.

txt_file : str
* Name of text file to save URLs to. Default is 'download_urls'.

save_to_txt : bool
* Boolean specifying if DEM URLs should be saved to a text file. Default is True.

download : bool
* Boolean specifying if DEM URLs retrieved should be downloaded immediately. Default is False.

#### Outputs
Text File
* Depending on passed parameters, may save URLs to a text file.

Folder
* Depending on passed parameters, may create a folder 'dem_tiles' for storing downloaded DEMs in.

DEM Files
* Depending on passed parameters, may download and store DEM files returned from fetch.

#### Notes
* If a shapefile and bounding box are both specified, the shapefile takes precedent
* GEOtiled currently only supports computation on GeoTIFF files, and it is not recommended to change the prod_format variable


### `get_extent(shape_file)`

#### Returns the bounding box (extents) of a shapefile.
This function extracts the extents of a shapefile. The extent is the upper left and lower right coordinates of the file.

#### Required Parameters
shape_file : str
* String specifying path to the shapefile.

#### Returns 
tuple of tuple
* Returns two tuples - the first being the upper left (x,y) coordinate, and the second being the lower right (x,y) coordinate


### `get_file_size(url)`

#### Return the size of a file, in bytes, retrieved from a URL.
This function uses the requests.head() function to read and return the 'Content-Length' header specified by a file found at a URL. This function is intended to be used with `download_files` to calculate download size before file downloads begin.

#### Required Parameters
url : str
* String representing URL where file is found.

#### Returns
int
* Size of the file specified at the URL in bytes. Returns 0 if file size can't be determined.


### `reproject(input_file, output_file, projection, cleanup=False)

#### Reprojects a GeoTIFF file to a specified projection.
This function reprojects a specified GeoTIFF file to a new projection changing the coordinate system representation, and the result is saved to a new file.

#### Required Parameters
input_file : str
* String specifying name of GeoTIFF file in data directory to reproject

output_file : str
* String specifying name of new reprojected GeoTIFF file to store in data directory

projection : str
* String specifying name of projection to use for reprojection. Can be a EPSG code (e.g. EPSG:4326) or the path to a WKT file.

#### Optional Parameters
cleanup : bool
* Boolean specifying if old file that was reprojected should be deleted after computation is complete. Default is False.


### `set_data_directory(path)`

#### Sets the path where data computed by GEOtiled will be stored.
This function sets the data directory where data will be searched for and generated in by GEOtiled functions.

#### Required Parameters
path : str
* A string that specifies working directory.

#### Notes
* It's good practice to run this function before executing anything else in the workflow.
* If not set, data will be searched for and stored in the working directory.


