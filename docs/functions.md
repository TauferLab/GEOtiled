# 🔧 Function Documentation

Here is where all useable functions from the GEOtiled library are documented.

## Functions

### `geotiled.build_stack(input_files, output_file, verbose=False)`  
Stacks multiple rasters into a single raster with multiple bands. The band order will be based off the input_files list order. 

**Parameters**  
input_files : List[str]  
&nbsp;&nbsp; Name/paths of rasters to stack together.  
output_file : str  
&nbsp;&nbsp; Name/path of file to store stacked raster data to.  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation of parameters (default is False).

**Notes**  
\- It is recommended that all rasters being stacked have a description for easier identification between bands in the stacked raster.

---

### `geotiled.crop_and_compute(input_file, parameter_list, tile_dimensions=None, num_tiles=None, compute_method="SAGA", convert_file=True, projection=5070, buffer_size=10, num_processes=None, cleanup=False, verbose=False)`  
Computes terrain parameters in parallel. This function handles folder creation, cropping of the elevation data into tiles, parallel computation of terrain parameters, and debuffering of computed terrain parameter tiles. Concurrent computation will, by default, use the max number of available cores for processing unless `num_processes` is specified.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path to input raster file with elevation data to compute terrain paramters from.  
parameter_list : List[str]  
&nbsp;&nbsp; List of terrain parameters to compute. If special keyword ‘all’ in list, will compute all terrain parameters for specific library.  
tile_dimensions : List[int], optional  
&nbsp;&nbsp; Column and row length of pixels to use for cropped tiles [x,y] (default is None).  
num_tiles : int, optional  
&nbsp;&nbsp; Will dynamically determine the size of the tiles if column_length and row_length is not specified (default is None).  
compute_method : str, optional  
&nbsp;&nbsp; API to use for computing terrain parameters (default is ‘SAGA’).  
convert_file : bool, optional  
&nbsp;&nbsp; Determine if files should be converted to SGRD when using SAGA as the compute method (default is True).  
projection : int, optional  
&nbsp;&nbsp; EPSG projection to set to metadata for created shapefiles (default is 5070).  
buffer_size : int, optional  
&nbsp;&nbsp; Number of buffer pixels to use for cropping (default is 10).  
num_processes : int, optional  
&nbsp;&nbsp; Number of concurrent processes to use for computing terrain parameters (default is None).  
cleanup: bool, optional  
&nbsp;&nbsp; Specifies if cropped files used for computing parameters should be deleted after computation (default is False).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation of parameters (default is False).

**Notes**  
\- It is recommended to use a tile size the produces an equivalent tile count that equals the number of concurrent processes being used for best performance results.  
\- It is important to set `projection` to the same projection of the `input_file`.

---

### `geotiled.crop_to_coordinates(input_file, output_file, upper_left, lower_right)`  
Crops a raster file based on a specified upper-left and lower-right coordinates. The upper_left and lower_right coordinates define the bounding box for cropping. Coordinates must be in the same projection as the raster.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path of the input raster file to crop.  
output_file : str  
&nbsp;&nbsp; Name of cropped raster to save the cropped data to.  
upper_left : tuple of float  
&nbsp;&nbsp; Float tuple specifying upper-left (x,y) coordinates to crop raster from.  
lower_right : tuple of float  
&nbsp;&nbsp; Float tuple specifying lower-right (x,y) coordinates to crop raster from.

---

### `geotiled.crop_to_region(input_file, output_file, codes)`  
Crop a raster file to the boundaries specified by multiple shapefiles. Shapefiles will be combined into a single geometry for cropping. The function automatically handles downloading shapefiles that are not present.

**Parameters**
input_file : str  
&nbsp;&nbsp; Name/path of input raster file to crop.  
output_file : str  
&nbsp;&nbsp; Name/path of output raster file cropped data will be written to.  
codes : List[str]  
&nbsp;&nbsp; List of shapefile codes that will outline the cropping region.

**Returns**  
shape_paths : List[str]  
&nbsp;&nbsp; Returns a list of paths to shapefiles used for cropping.

---

### `geotiled.crop_to_valid_data(input_file, output_file, projection="EPSG:4269", block_size=512)`  
Crops a raster file to the extents of valid data, where valid data is a row or column of data that contains at least one non-nan value (i.e., it crops borders of the raster file containing rows or columns of only nan values). Blocking is used to help minimize RAM usage, and should be adjusted accordingly.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path of raster file to crop.  
output_file : str  
&nbsp;&nbsp; Name/path of raster file to save cropped data to.  
projection : str, optional  
&nbsp;&nbsp; Name of projection to reference for translation of new file (default is 'EPSG:4269').  
block_size : int, optional  
&nbsp;&nbsp; Block size to use when computing extents (default is 512).

---

### `geotiled.crop_to_window(input_file, output_file, window)`  
Crops a raster file to a specific window where the window references 2D matrix indices.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path to input file to be cropped.  
output_file : str  
&nbsp;&nbsp; Name/path to file to write the cropped data.  
window : List[int], Tuple(int)  
&nbsp;&nbsp; List or tuple of format [left_x, top_y, width, height] where left_x and top_y are matrix indices of the upper-left corner of the cropping window, and width and height specify the dimensions of the cropping window in pixels.

---

### `geotiled.download_files(download_list, download_folder)`  
Download file(s) from a list of URLs simultaneously using threading and showcases download progress via a tqdm progress bar.

**Parameters**  
download_list : str, List[str]  
&nbsp;&nbsp; The name/path of a text file. This file should contain download URLs separated by newlines. A list of download URLs.  
download_folder : str  
&nbsp;&nbsp; Name/path of the folder where the downloaded files will be stored.

---

### `geotiled.download_shapefiles(codes)`  
Downloads specified shapefile(s) off the USGS webpage and stores them in a shapefiles folder located in the working directory. All codes passed should be valid US state abbreviations. It will skip downloading already existing shapefiles.

**Parameters**  
codes : str, List[str]  
&nbsp;&nbsp; Comma-separated string of valid codes. List of valid codes.

---

### `geotiled.extract_raster(csv_file, raster_file, band_names, verbose=False)`  
Uploads raster data to a CSV that already has x,y coordinates. Only data from the raster that has matching x,y coordinates are written to the CSV. Order of band_names should correlate to order of the raster’s bands. If no x,y coordinates in the CSV correlate to those of the raster file, incorrect or no values will be extracted.

**Parameters**  
csv_file : str  
&nbsp;&nbsp; Name/path to CSV file to read/write to.  
raster_file : str  
&nbsp;&nbsp; Name/path to raster file to read to values from.  
band_names : List[str], optional  
&nbsp;&nbsp; Names of columns correlating to raster bands. Order matters.  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation of parameters (default is False).

---

### `geotiled.fetch_dems(shapefile=None, bbox={"xmax": -83.815, "xmin": -84.0387, "ymax": 36.04, "ymin": 35.86}, dataset="30m", txt_file="download_urls.txt", save_to_txt=True, download_folder="dem_tiles", download=False, verbose=False)`  
Queries USGS National Map API to fetch DEM data URLs using specified filters and can either save the list of URLs to a text file and/or download from the list of URLs immediately.

**Parameters**  
shapefile : str, optional  
&nbsp;&nbsp; Code of shapefile with which a bounding box will be generated (default is None). Overrides the bbox parameter if set.  
bbox : dict, optional  
&nbsp;&nbsp; Bounding box coordinates to query for DEM data (default is {'xmin': -84.0387, 'ymin': 35.86, 'xmax': -83.815, 'ymax': 36.04}).  
dataset : str, optional  
&nbsp;&nbsp; Resolution of DEM data to download (default is '30m').  
txt_file : str, optional  
&nbsp;&nbsp; Name of text file to save URLs to (default is 'download_urls.txt').  
save_to_txt : bool, optional  
&nbsp;&nbsp; Allow DEM URLs to be saved to a text file (default is True).  
download_folder : str, optional  
&nbsp;&nbsp; Name of the download folder to store downloaded DEMs in (default is 'dem_tiles').  
download : bool, optional  
&nbsp;&nbsp; Allow DEM URLs retrieved to be immediately downloaded (default is False).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track download (default is False).

**Notes**  
\- Available resolutions to download can be found in the [Data Guide](./data.md).

---

### `geotiled.merge_shapefiles(input_folder, output_file, cleanup=False, verbose=False)`  
This function merges multiple shapefiles together into a single shapefile. Shapefiles provided should be in the .shp format.

**Parameters**  
input_folder : str  
&nbsp;&nbsp; Name of folder where shapefiles to merge are stored.  
output_file : str  
&nbsp;&nbsp; Name of output file that has merged shapefiles.  
cleanup : bool, optional  
&nbsp;&nbsp; Determine if files from input folder should be deleted after computation (default is False).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation (default is False).

**Notes**  
\- Ensure that all shapefiles you wish to merge are located in the `input_folder`.

---

### `geotiled.mosaic_rasters(input_folder, output_file, description=None, cleanup=False, verbose=False)`  
Mosaics together raster files into a single raster. The raster will only contain one band, and the user has the option to update the description of the band.

**Parameters**  
input_folder : str  
&nbsp;&nbsp; Name/path to the folder containing rasters to mosaic.  
output_file : str  
&nbsp;&nbsp; Name/path to the file to write the mosaicked rasters data to.  
description : str, optional  
&nbsp;&nbsp; Description to add to mosaicked raster’s band (default is None).  
cleanup : bool, optional  
&nbsp;&nbsp; Determines if files from input_folder should be deleted after mosaic is complete (default is False).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation (default is False).

**Notes**  
\- Ensure that all rasters you wish to merge are located in the `input_folder`.

---

### `geotiled.plot_raster(input_file, plot_title=None, reproject_gcs=False, projection='EPSG:4269', shapefiles=None, remove_nans=False, downsample=1, dpi=150, cmap='inferno', nancolor='white', ztype='Z', zunit=None, xyunit=None, vmin=None, vmax=None, bordercolor='black', borderlinewidth=1.5, clean=False, save_to_img=None, verbose=False)`  
Plots a raster with user-specified modifications. The user also has the option to save the plot to a file.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path of the raster file to plot.  
plot_title : str  
&nbsp;&nbsp; Title for plot (default is None).  
reproject_gcs : bool, optional  
&nbsp;&nbsp; Determine if raster should be reprojected before visualization (Default is False).  
projection : str, optional  
&nbsp;&nbsp; Projection to change raster to if reprojection is specified (default is ‘EPSG:4269’).  
shapefiles : List[str]  
&nbsp;&nbsp; List of valid shapefile codes to crop raster data to before visualization (default is None).  
remove_nans : bool, optional  
&nbsp;&nbsp; Determine if nan-borders on data should be removed before plotting (default is False).  
downsample : int, optional  
&nbsp;&nbsp; Factor to downsample the raster by (default is 1).  
dpi : int, optional  
&nbsp;&nbsp; Resolution in dots per inch for the figure. Default is 150.  
cmap : str, optional  
&nbsp;&nbsp; Colormap used for visualization (default is ‘inferno’).  
nancolor : str, optional  
&nbsp;&nbsp; Color to use for NaN values (default is ‘white’).  
ztype : str, optional  
&nbsp;&nbsp; Data that is represented by the z-axis (default is ‘Z’).  
zunit : str, optional  
&nbsp;&nbsp; Units for the data values (z-axis) which is inferred unless specified (default is None).  
xyunit : str, optional  
&nbsp;&nbsp; Units for the x and y axes which is inferred unless specified (default is None).  
vmin : float, optional  
&nbsp;&nbsp; Value of lower bound for coloring on plot which is the minimum of the data unless specified (default is None).  
vmax : float, optional  
&nbsp;&nbsp; Value of upper bound for coloring on plot which is the maximum of the data unless specified (default is None).  
bordercolor : str, optional  
&nbsp;&nbsp; Color for the shapefile boundary if one is used for cropping (default is ‘black’).  
borderlinewidth : float, optional  
&nbsp;&nbsp; Line width for the shapefile boundary if one is used for cropping (default is 1.5).  
clean : bool, optional  
&nbsp;&nbsp; Determine whether to plot only the image (default is False).  
save_to_img : str, optional  
&nbsp;&nbsp; Name/path to image file to save plot to (default is None).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation of parameters (default is False).  

**Outputs**  
Image  
&nbsp;&nbsp; Displays the plot of a raster with specified modifications.

**Returns**  
raster_array : np.ndarray  
&nbsp;&nbsp; The data, in array format, of the plotted raster.

**Notes**  
\- Downsampling is recommended for larger raster files to ensure efficient plotting.

---

### `geotiled.plot_shapefile(input_file, plot_title='', reproject_gcs=False, projection=4269, crop_to_shape=None, save_to_img=None, verbose=False)`  
Visualizes a shapefile in a plot. The projection can be changed, and the data can be cropped to the bounds of another shapefile. The user also has the option to save the plot to a file.

**Parameters**
input_file : str  
&nbsp;&nbsp; Name/path to the shapefile to plot.  
plot_title : str, optional  
&nbsp;&nbsp; Title to assign to plot (default is ‘’).  
reproject_gcs : bool, optional  
&nbsp;&nbsp; Determine whether to reproject data before visualization (default is False).  
projection : int, optional  
&nbsp;&nbsp; If reprojecting data, specifies which EPSG to reproject to (default is 4269).  
crop_to_shape : str, optional  
&nbsp;&nbsp; State shapefile to crop data to (default is None).  
save_to_img : str, optional  
&nbsp;&nbsp; Name/path of file to save plot to (default is None).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation of parameters (default is False).

**Outputs**  
Image  
&nbsp;&nbsp; Creates a plot for the shapefile.

---

### `geotiled.print_computable_parameters()`  
Outputs all computable terrain parameters with GEOtiled. A special indicator ‘(G)’ indicates it is computable with GDAL.

---

### `geotiled.print_region_codes()`  
Outputs state abbreviations and their correlating states. 

---

### `geotiled.reproject(input_file, output_file, projection, cleanup=False, verbose=False)`  
Reprojects a raster file to a specified projection where the result is saved to a new file. Multithreading is utilized to improve performance.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path of raster to reproject.  
output_file : str  
&nbsp;&nbsp; Name/path of file to write reprojected data to.  
projection : str  
&nbsp;&nbsp; Projection to use for reprojection. Can be a EPSG code (e.g. EPSG:4269) or the name/path to a WKT file.  
cleanup : bool, optional  
&nbsp;&nbsp; Determine if input_file should be deleted after reprojection is complete (default is False).  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation (default is False).

**Notes**  
\- It is important to use this for preprocessing DEM data to have units of meters.

---

### `geotiled.set_working_directory(path)`  
Creates needed directories then sets the specified path as the working directory. It is important to run this before anything else to ensure data is stored in a desired location.

**Parameters**  
path : str  
&nbsp;&nbsp; Working directory to store data in.

---

### `geotiled.tif2csv(input_file, output_file, band_names, verbose=False)`  
Converts a GeoTIFF file and converts it into the CSV format. The CSV columns are the x coordinate, y coordinate, and raster values. NaN values indicate no data found at a particular index. If the GeoTIFF has multiple bands, the user can set the names of subsequent columns in the CSV containing each bands value for each x,y coordinate using the band_names parameter.

**Parameters**  
input_file : str  
&nbsp;&nbsp; Name/path to GeoTIFF file to convert to a CSV.  
output_file : str  
&nbsp;&nbsp; Name/path to CSV file to save CSV data to.  
band_names : List[str], optional  
&nbsp;&nbsp; Names of columns correlating to raster bands. Order matters.  
verbose : bool, optional  
&nbsp;&nbsp; Determine if additional print statements should be used to track computation of parameters (default is False).