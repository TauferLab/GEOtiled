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


### `build_mosaic_filtered(input_folder, output_file, cleanup=False)`

#### Builds mosaic from multiple GeoTIFF files that were cropped with buffer regions.
This function is similar to the `build_mosaic` function but handles mosaicking together GeoTIFF files that were split to includes buffer regions and averages points in the buffer regions together when merging.

#### Required Parameters
input_folder : str
* String specifying name of folder in data directory where GeoTIFF files to mosaic together are located.

output_file : str
* String specifying name of mosaicked file produced.

#### Optional Parameters
cleanup : bool
* Boolean specifying if tiles in input_folder used for computation should be deleted after computation is complete. Default is False.

#### Outputs
GeoTIFF File
* Produces a GeoTIFF file that is the mosaic of all GeoTIFF files from the specified input folder.


### `build_stack(input_files, output_file)`

#### Stacks multiple GeoTIFF files into a single GeoTIFF with multiple bands.
This function takes multiple GeoTIFF files and combines them into one file represented by different bands. This is useful when multiple datasets need to be represented as one.

#### Required Parameters
input_files : str list
* String list specifying paths to GeoTIFF files to be stacked together.

output_file : str
* String specifying path to store stacked raster file.

#### Outputs
GeoTIFF File
* Generates a single GeoTIFF file with multiple bands with the path it is stored in being establed by the variable 'output_file'.

#### Notes
* The band order will be based off the 'input_files' list order


### `change_raster_format(input_file, output_file, raster_format)`

#### Convert format of a specified raster file.
This function uses GDAL functions to convert the file type of specified raster data.

#### Required Parameters
input_file : str
* String specifying the path to the input file to be converted.

output_file : str
* String specifying the path to the output file that is converted.

raster_format : str
* GDAL supported string specifying the raster format to convert the input file to.

#### Outputs
Raster File
* Generates a raster of the format specified by 'raster_format'.

#### Notes
* Supported formats can be found of the [GDAL website](https://gdal.org/drivers/raster/index.html)


### `compute_geotiled(input_folder, param_list, num_procs, cleanup=False)`

#### Configures the multiprocessing pool for GEOtiled to begin computing terrain parameters.
This function utilizes the multiprocessing library to allow for computation of parameters on different elevation GeoTIFF files at the same time.

#### Required Parameters
input_folder : str
* String containing the name of the folder in the data directory containing DEM elevation GeoTIFF files.

param_list : str list
* String list containing codes for terrain parameters to compute. The 'all' keyword will compute all params.

num_procs : int
* Integer specifying the number of python instances to use for multiprocessing.

#### Optional Parameters
cleanup : bool
* Boolean specifying if elevation files used to compute parameters should be deleted after computation. Default is False.

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


### `crop_coord(input_file, output_file, upper_left, lower_right)`

#### Crops raster file to a specific region given specified coordinates.
This function uses GDAL functions to crop a raster file based on a specified upper-left and lower-right coordinate.

#### Required Parameters
input_file : str
* String specifying path to input raster to be cropped.

output_file : str
* String specifying path to store cropped raster to.

upper_left : float tuple
* Float tuple specifying upper-left (x,y) coordinates to crop raster from.

lower_right : float tuple
* Float tuple specifying lower-right (x,y) coordinates to crop raster from.

#### Outputs
Cropped GeoTIFF File
* Will generate the cropped GeoTIFF file at the location specified by 'output_file'. 

#### Error States
* An error will be thrown if coordinates specified fall outside the range of the input raster's bounds.


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


### `crop_region(input_file, shape_file, output_file)`

#### Crops a raster file based off the region defined by a shapefile.
This function uses a GDAL function to crop a raster file according to boundaries specified by a shapefile.

#### Required Parameters
input_file : str
* String specifying path to GeoTIFF file to be cropped.

shape_file : str
* String specifying path to shapefile to use for cropping.

output_file : str
* String specifying path where cropped raster will be saved.

#### Outputs
Cropped GeoTIFF File
* A GeoTIFF File with data cropped to the bounds of the shapefile is created.

#### Error States
* Error will generate if bounds of shapefile exceed input raster's limits


### `crop_to_valid_data(input_file, output_file, block_size=512)`

#### Crops a border region of NaN values from a GeoTIFF file.
This function uses blocking to scan through a GeoTIFF file to determine the extent of valid data and crops excess Nan values at the borders of the spatial data.

#### Required Parameters
input_file : str
* String specifying path to GeoTIFF file to crop.

output_file : str
* String specifying path to GeoTIFF file to save cropped data to.

#### Optional Parameters
block_size : int
* Integer specifying the block size to use when computing extents. Default is 512.

#### Outputs
GeoTIFF file
* Output is a cropped GeoTIFF file with the path specified by 'output_file'

#### Notes
* Blocking helps minimize RAM usage, and should be adjusted as needed to help improve performance.


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


### `extract_raster(csv_file, raster_file, band_names)`

#### Extracts raster values and stores them with correlating coordinates found in a CSV file.
This function reads raster values from a file and stores the value with correlating x and y coordinates in a CSV file.

#### Required Parameters
csv_file : str
* String specifying path to CSV file for input.

raster_file : str
* String specifying path to raster file to read raster values from.

band_names : str list
* String list specifying names of bands to extract values from in input raster.

#### Outputs
CSV File
* An updated CSV file with raster values correlating to x and y coordinates will be produced.

#### Notes
* CSV file must already be appropiately formatted with x and y coordinates.
* Order of 'band_names' should correlate to order of band names in specified raster file.
* If coordinates in CSV file are outside of bounds of raster, incorrect or no values will be extracted.


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


### `generate_img(tif, cmap='inferno', dpi=150, downsample=1, verbose=False, clean=False, title=None, nancolor='green', ztype="Z", zunit=None, xyunit=None, vmin=None, vmax=None, reproject_gcs=False, shp_files=None, crop_shp=False, bordercolor="black", borderlinewidth=1.5, saveDir = None)`

#### Plots a GeoTIFF image using matplotlib.
This function is meant to visualize GeoTIFF with optional parameters to ensure data visualization is made easily discernable and customizable.

#### Required Parameters
tif : str
* String specifying name of the GeoTIFF file in the data directory to plot.

#### Optional Parameters
cmap : str
* Colormap used for visualization. Default is 'inferno'.
  
dpi : int
* Resolution in dots per inch for the figure. Default is 150.
    
downsample : int
* Factor to downsample the image by. Default is 10.
    
verbose : bool
* If True, print geotransform and spatial reference details. Default is False.
    
clean : bool
* If True, no extra data will be shown besides the plot image. Default is False.
    
title : str
* Title for the plot. Default will display the projection name.
    
nancolor : str
* Color to use for NaN values. Default is 'green'.
    
ztype : str
* Data that is represented by the z-axis. Default is 'Z'.
    
zunit : str
* Units for the data values (z-axis). Default is None and inferred from spatial reference.
    
xyunit : str
* Units for the x and y axes. Default is None and inferred from spatial reference.
    
vmin : float
* Value of lower bound for coloring on plot. Will be whatever the min value of the data is if none is specified.
    
vmax : float
* Value of upper bound for coloring on plot. Will be whatever the max value of the data is if none is specified.
    
reproject_gcs : bool
* Reproject a given raster from a projected coordinate system (PCS) into a geographic coordinate system (GCS).
    
shp_files : str list
* Comma-seperated list of strings with shape file codes to use for cropping. Default is None.
    
crop_shp : bool
* Flag to indicate if the shapefiles should be used for cropping. Default is False.
    
bordercolor : str
* Color for the shapefile boundary. Default is "black".
    
borderlinewidth : float
* Line width for the shapefile boundary. Default is 1.5.

#### Outputs
Image
* Displays image visualizing inputed GeoTIFF data with specified parameters.

#### Returns
raster_array : np.ndarray
* Returns the raster array that was used for visualization. 

#### Notes
* Alternative colormaps can be found in the [matplotlib documentation](https://matplotlib.org/stable/users/explain/colors/colormaps.html).
* Graph currently only tested for visualization with Jupyter Notebooks. 
* Using 'shp_file' without setting 'crop_shp' will allow you to plot the outline of the shapefile without actually cropping anything.


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


### `reproject(input_file, output_file, projection, cleanup=False)`

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


### `tif2csv(raster_file, band_names=['elevation'], output_file='params.csv')`

#### Converts a raster file into a CSV file.
This function reads values from a GeoTIFF file and converts it into the CSV format. The CSV columns are the x coordinate, y coordinate, and raster value.

#### Required Parameters
raster_file : str
* String specifying path to raster file to be stored as a CSV.

#### Optional Parameters
band_names : str list
* String list specifying names of GeoTIFF bands to pull data from.

output_file : str
* String specifying path to CSV storing converted raster values.

#### Outputs
CSV File
* Generates CSV file with converted data from raster file.

#### Notes
* NaN values indicate no data found at a particular coordinate
* Only band names provided will be read into the CSV. If missing band names from GeoTIFF file, CSV may contain columns without headers or be missing data.
