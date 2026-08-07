# 🔧 Function Documentation

Here is where all useable functions from the GEOtiled library are documented.

## 🌎 North American Ecoregions

### `ECOREGION_TERRAIN_PARAMETERS`

```python
geotiled.ECOREGION_TERRAIN_PARAMETERS
```

The standard ecoregion output bundle contains 16 SAGA terrain parameters. It produces 14 raster parameters and the `channel_network` and `drainage_basins` vector parameters. Cropping the source DEM as `elevation.tif` creates 17 final deliverables for an ecoregion.

---

### `build_stack()`

```python
geotiled.build_stack(input_files, output_file, verbose=False)
```

Stacks multiple rasters into a single raster with multiple bands. The band order will be based off the input_files list order.

**Parameters**

| Name        | Type        | Description                                                     | Default |
| ----------- | ----------- | --------------------------------------------------------------- | ------- |
| input_files | `List[str]` | Names/paths of rasters to stack together.                       | —       |
| output_file | `str`       | Name/path of file to store stacked raster data.                 | —       |
| verbose     | `bool`      | Whether to print additional progress messages during execution. | `False` |

> **_NOTE:_** It is recommended that all rasters being stacked have a description for easier identification between bands in the stacked raster.

---

### `crop_and_compute()`

```python
geotiled.crop_and_compute(input_file, parameter_list, tile_dimensions=None, num_tiles=None, compute_method="SAGA", convert_file=True, projection=5070, buffer_size=10, num_processes=None, cleanup=False, verbose=False)
```

Computes terrain parameters in parallel. This function handles folder creation, cropping of the elevation data into tiles, parallel computation of terrain parameters, and debuffering of computed terrain parameter tiles. Concurrent computation will, by default, use the max number of available cores for processing unless `num_processes` is specified.

**Parameters**

| Name            | Type        | Description                                                                                                                        | Default  |
| --------------- | ----------- | ---------------------------------------------------------------------------------------------------------------------------------- | -------- |
| input_file      | `str`       | Name/path to input raster file with elevation data to compute terrain paramters from.                                              | —        |
| parameter_list  | `List[str]` | List of terrain parameters to compute. If special keyword ‘all’ in list, will compute all terrain parameters for specific library. | —        |
| tile_dimensions | `List[int]` | Column and row length of pixels to use for cropped tiles [x,y].                                                                    | `None`   |
| num_tiles       | `int`       | Will dynamically determine the size of the tiles if column_length and row_length is not specified.                                 | `None`   |
| compute_method  | `str`       | API to use for computing terrain parameters.                                                                                       | `'SAGA'` |
| convert_file    | `bool`      | Determine if files should be converted to SGRD when using SAGA as the compute method.                                              | `True`   |
| projection      | `int`       | EPSG projection to set to metadata for created shapefiles.                                                                         | `5070`   |
| buffer_size     | `int`       | Number of buffer pixels to use for cropping.                                                                                       | `10`     |
| num_processes   | `int`       | Number of concurrent processes to use for computing terrain parameters.                                                            | `None`   |
| cleanup         | `bool`      | Specifies if cropped files used for computing parameters should be deleted after computation.                                      | `False`  |
| verbose         | `bool`      | Determine if additional print statements should be used to track computation of parameters.                                        | `False`  |

> **NOTE:** It is recommended to use a tile size the produces an equivalent tile count that equals the number of concurrent processes being used for best performance results.  
> **NOTE:** It is important to set `projection` to the same projection of the `input_file`.

---

### `crop_to_coordinates()`

```python
geotiled.crop_to_coordinates(input_file, output_file, upper_left, lower_right)
```

Crops a raster file based on a specified upper-left and lower-right coordinates. The upper_left and lower_right coordinates define the bounding box for cropping. Coordinates must be in the same projection as the raster.

**Parameters**

| Name        | Type             | Description                                                               | Default |
| ----------- | ---------------- | ------------------------------------------------------------------------- | ------- |
| input_file  | `str`            | Name/path of the input raster file to crop.                               | —       |
| output_file | `str`            | Name of cropped raster to save the cropped data to.                       | —       |
| upper_left  | `tuple of float` | Float tuple specifying upper-left (x,y) coordinates to crop raster from.  | —       |
| lower_right | `tuple of float` | Float tuple specifying lower-right (x,y) coordinates to crop raster from. | —       |

---

### `crop_to_ecoregion()`

```python
geotiled.crop_to_ecoregion(input_file, output_file, ecoregion)
```

Crops a raster to a dissolved EPA Level I, II, or III North American ecoregion boundary. The level is inferred from the string code, and the corresponding EPA archive is downloaded and cached beneath the working directory's `shapefiles/ecoregions` folder when needed. GDAL transforms the cutline to the raster's coordinate reference system during the crop.

**Parameters**

| Name        | Type  | Description                                               | Default |
| ----------- | ----- | --------------------------------------------------------- | ------- |
| input_file  | `str` | Name/path of the raster to crop.                          | —       |
| output_file | `str` | Name/path where the cropped raster will be written.       | —       |
| ecoregion   | `str` | EPA Level I–III code, such as `11`, `11.1`, or `11.1.1`. | —       |

**Returns**

| Type  | Description                      |
| ----- | -------------------------------- |
| `str` | The supplied output raster path. |

**Raises**

| Type           | Description                                                              |
| -------------- | ------------------------------------------------------------------------ |
| `ValueError`   | The code is malformed, unsupported, `0`/`WATER`, or absent from EPA data. |
| `RuntimeError` | The EPA boundary cannot be loaded or GDAL cannot create the crop.        |

---

### `crop_to_region()`

```python
geotiled.crop_to_region(input_file, output_file, codes)
```

Crop a raster file to the boundaries specified by multiple shapefiles. Shapefiles will be combined into a single geometry for cropping. The function automatically handles downloading shapefiles that are not present.

**Parameters**

| Name        | Type        | Description                                                      | Default |
| ----------- | ----------- | ---------------------------------------------------------------- | ------- |
| input_file  | `str`       | Name/path of input raster file to crop.                          | —       |
| output_file | `str`       | Name/path of output raster file cropped data will be written to. | —       |
| codes       | `List[str]` | List of shapefile codes that will outline the cropping region.   | —       |

**Returns**

| Name        | Type        | Description                                              | Default |
| ----------- | ----------- | -------------------------------------------------------- | ------- |
| shape_paths | `List[str]` | Returns a list of paths to shapefiles used for cropping. | —       |

---

### `crop_to_valid_data()`

```python
geotiled.crop_to_valid_data(input_file, output_file, projection="EPSG:4269", block_size=512)
```

Crops a raster file to the extents of valid data, where valid data is a row or column of data that contains at least one non-nan value (i.e., it crops borders of the raster file containing rows or columns of only nan values). Blocking is used to help minimize RAM usage, and should be adjusted accordingly.

**Parameters**

| Name        | Type  | Description                                                  | Default       |
| ----------- | ----- | ------------------------------------------------------------ | ------------- |
| input_file  | `str` | Name/path of raster file to crop.                            | —             |
| output_file | `str` | Name/path of raster file to save cropped data to.            | —             |
| projection  | `str` | Name of projection to reference for translation of new file. | `'EPSG:4269'` |
| block_size  | `int` | Block size to use when computing extents.                    | `512`         |

---

### `crop_to_window()`

```python
geotiled.crop_to_window(input_file, output_file, window)
```

Crops a raster file to a specific window where the window references 2D matrix indices.

**Parameters**

| Name        | Type               | Description                                                                                                                                                                                                                     | Default |
| ----------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- |
| input_file  | `str`              | Name/path to input file to be cropped.                                                                                                                                                                                          | —       |
| output_file | `str`              | Name/path to file to write the cropped data.                                                                                                                                                                                    | —       |
| window      | `List[int], Tuple` | List or tuple of format [left_x, top_y, width, height] where left_x and top_y are matrix indices of the upper-left corner of the cropping window, and width and height specify the dimensions of the cropping window in pixels. | —       |

---

### `download_files()`

```python
geotiled.download_files(download_list, download_folder)
```

Download file(s) from a list of URLs simultaneously using threading and showcases download progress via a tqdm progress bar.

**Parameters**

| Name            | Type               | Description                                                                                                          | Default |
| --------------- | ------------------ | -------------------------------------------------------------------------------------------------------------------- | ------- |
| download_list   | `str`, `List[str]` | The name/path of a text file. This file should contain download URLs separated by newlines. A list of download URLs. | —       |
| download_folder | `str`              | Name/path of the folder where the downloaded files will be stored.                                                   | —       |

---

### `download_shapefiles()`

```python
geotiled.download_shapefiles(codes)
```

Downloads specified shapefile(s) off the USGS webpage and stores them in a shapefiles folder located in the working directory. All codes passed should be valid US state abbreviations. It will skip downloading already existing shapefiles.

**Parameters**

| Name  | Type               | Description                                                 | Default |
| ----- | ------------------ | ----------------------------------------------------------- | ------- |
| codes | `str`, `List[str]` | Comma-separated string of valid codes. List of valid codes. | —       |

---

### `extract_raster()`

```python
geotiled.extract_raster(csv_file, raster_file, band_names, verbose=False)
```

Uploads raster data to a CSV that already has x,y coordinates. Only data from the raster that has matching x,y coordinates are written to the CSV. Order of band_names should correlate to order of the raster’s bands. If no x,y coordinates in the CSV correlate to those of the raster file, incorrect or no values will be extracted.

**Parameters**

| Name        | Type        | Description                                                                                 | Default |
| ----------- | ----------- | ------------------------------------------------------------------------------------------- | ------- |
| csv_file    | `str`       | Name/path to CSV file to read/write to.                                                     | —       |
| raster_file | `str`       | Name/path to raster file to read to values from.                                            | —       |
| band_names  | `List[str]` | Names of columns correlating to raster bands. Order matters.                                | —       |
| verbose     | `bool`      | Determine if additional print statements should be used to track computation of parameters. | `False` |

---

### `fetch_dems()`

```python
geotiled.fetch_dems(shapefile=None, bbox={"xmin": -84.0387, "ymin": 35.86, "xmax": -83.815, "ymax": 36.04}, dataset="30m", txt_file="download_urls.txt", save_to_txt=True, download_folder="dem_tiles", download=False, verbose=False, ecoregion=None)
```

Queries the USGS National Map API for DEM URLs and optionally saves or downloads the matching files. Supply a state `shapefile`, an explicit `bbox`, or an EPA North American ecoregion code.

**Parameters**

| Name            | Type   | Description                                                                                                                    | Default                                                             |
| --------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------- |
| shapefile       | `str`  | State code used to generate the query bounds. Cannot be combined with `ecoregion`.                                            | `None`                                                              |
| bbox            | `dict` | Bounding box coordinates to query when neither `shapefile` nor `ecoregion` is supplied.                                     | `{'xmin': -84.0387, 'ymin': 35.86, 'xmax': -83.815, 'ymax': 36.04}` |
| dataset         | `str`  | Resolution code for DEM data. See the [Data Guide](./data.md).                                                                 | `'30m'`                                                             |
| txt_file        | `str`  | Name/path of the text file where URLs are saved.                                                                               | `'download_urls.txt'`                                               |
| save_to_txt     | `bool` | Whether to save the retrieved URLs to `txt_file`.                                                                             | `True`                                                              |
| download_folder | `str`  | Name/path of the folder where downloaded DEMs are stored.                                                                      | `'dem_tiles'`                                                       |
| download        | `bool` | Whether to download the retrieved DEMs immediately.                                                                            | `False`                                                             |
| verbose         | `bool` | Whether to print additional download progress.                                                                                 | `False`                                                             |
| ecoregion       | `str`  | EPA Level I–III code. The number of dot-separated components determines the level.                                             | `None`                                                              |

For ecoregions, the matching EPA archive is cached beneath the working directory's `shapefiles/ecoregions` folder. Broad boundaries are divided into intersecting 5-degree query windows. Each USGS request is paginated, and duplicate download URLs are removed.

**Raises**

| Type         | Description                                                               |
| ------------ | ------------------------------------------------------------------------- |
| `ValueError` | The ecoregion code is invalid or is supplied together with `shapefile`. |

---

### `merge_shapefiles()`

```python
geotiled.merge_shapefiles(input_folder, output_file, cleanup=False, verbose=False, ecoregion=None)
```

Merges multiple shapefiles into one `.shp` file or ZIP archive. When an ecoregion is provided, the merged features are reprojected as needed and clipped to its dissolved EPA boundary. A `.zip` output contains the `.shp`, `.shx`, `.dbf`, `.prj`, and optional `.cpg` components.

**Parameters**

| Name         | Type   | Description                                                                   | Default |
| ------------ | ------ | ----------------------------------------------------------------------------- | ------- |
| input_folder | `str`  | Name/path of the folder containing shapefiles to merge.                       | —       |
| output_file  | `str`  | Name/path of the output `.shp` or `.zip` file.                               | —       |
| cleanup      | `bool` | Whether to delete `input_folder` after a successful merge.                    | `False` |
| verbose      | `bool` | Whether to print additional progress messages.                                | `False` |
| ecoregion    | `str`  | EPA Level I–III code used to clip the merged features.                        | `None`  |

**Raises**

| Type         | Description                                                  |
| ------------ | ------------------------------------------------------------ |
| `ValueError` | The ecoregion code is invalid or absent from the EPA layer. |

---

### `mosaic_rasters()`

```python
geotiled.mosaic_rasters(input_folder, output_file, description=None, cleanup=False, verbose=False)
```

Mosaics together raster files into a single raster. The raster will only contain one band, and the user has the option to update the description of the band.

**Parameters**

| Name         | Type   | Description                                                                       | Default |
| ------------ | ------ | --------------------------------------------------------------------------------- | ------- |
| input_folder | `str`  | Name/path to the folder containing rasters to mosaic.                             | —       |
| output_file  | `str`  | Name/path to the file to write the mosaicked rasters data to.                     | —       |
| description  | `str`  | Description to add to mosaicked raster’s band.                                    | `None`  |
| cleanup      | `bool` | Determines if files from input_folder should be deleted after mosaic is complete. | `False` |
| verbose      | `bool` | Determine if additional print statements should be used to track computation.     | `False` |

> **NOTE:** Ensure that all rasters you wish to merge are located in the `input_folder`.

---

### `plot_raster()`

```python
geotiled.plot_raster(input_file, plot_title=None, reproject_gcs=False, projection='EPSG:4269', shapefiles=None, remove_nans=False, downsample=1, dpi=150, cmap='inferno', nancolor='white', ztype='Z', zunit=None, xyunit=None, vmin=None, vmax=None, bordercolor='black', borderlinewidth=1.5, clean=False, save_to_img=None, verbose=False)
```

Plots a raster with user-specified modifications. The user also has the option to save the plot to a file.

**Parameters**

| Name            | Type        | Description                                                                                  | Default       |
| --------------- | ----------- | -------------------------------------------------------------------------------------------- | ------------- |
| input_file      | `str`       | Name/path of the raster file to plot.                                                        | —             |
| plot_title      | `str`       | Title for plot.                                                                              | `None`        |
| reproject_gcs   | `bool`      | Determine if raster should be reprojected before visualization.                              | `False`       |
| projection      | `str`       | Projection to change raster to if reprojection is specified.                                 | `'EPSG:4269'` |
| shapefiles      | `List[str]` | List of valid shapefile codes to crop raster data to before visualization.                   | `None`        |
| remove_nans     | `bool`      | Determine if nan-borders on data should be removed before plotting.                          | `False`       |
| downsample      | `int`       | Factor to downsample the raster by.                                                          | `1`           |
| dpi             | `int`       | Resolution in dots per inch for the figure.                                                  | `150`         |
| cmap            | `str`       | Colormap used for visualization.                                                             | `'inferno'`   |
| nancolor        | `str`       | Color to use for NaN values.                                                                 | `'white'`     |
| ztype           | `str`       | Data that is represented by the z-axis.                                                      | `'Z'`         |
| zunit           | `str`       | Units for the data values (z-axis) which is inferred unless specified.                       | `None`        |
| xyunit          | `str`       | Units for the x and y axes which is inferred unless specified.                               | `None`        |
| vmin            | `float`     | Value of lower bound for coloring on plot which is the minimum of the data unless specified. | `None`        |
| vmax            | `float`     | Value of upper bound for coloring on plot which is the maximum of the data unless specified. | `None`        |
| bordercolor     | `str`       | Color for the shapefile boundary if one is used for cropping.                                | `'black'`     |
| borderlinewidth | `float`     | Line width for the shapefile boundary if one is used for cropping.                           | `1.5`         |
| clean           | `bool`      | Determine whether to plot only the image.                                                    | `False`       |
| save_to_img     | `str`       | Name/path to image file to save plot to.                                                     | `None`        |
| verbose         | `bool`      | Determine if additional print statements should be used to track computation of parameters.  | `False`       |

**Returns**

| Name         | Type         | Description                                       | Default |
| ------------ | ------------ | ------------------------------------------------- | ------- |
| raster_array | `np.ndarray` | The data, in array format, of the plotted raster. | —       |

**Outputs**

| Name  | Type | Description                                                 | Default |
| ----- | ---- | ----------------------------------------------------------- | ------- |
| Image | —    | Displays the plot of a raster with specified modifications. | —       |

> **NOTE:** Downsampling is recommended for larger raster files to ensure efficient plotting.

---

### `plot_shapefile()`

```python
geotiled.plot_shapefile(input_file, plot_title='', reproject_gcs=False, projection=4269, crop_to_shape=None, save_to_img=None, verbose=False)
```

Visualizes a shapefile in a plot. The projection can be changed, and the data can be cropped to the bounds of another shapefile. The user also has the option to save the plot to a file.

**Parameters**

| Name          | Type   | Description                                                                                 | Default |
| ------------- | ------ | ------------------------------------------------------------------------------------------- | ------- |
| input_file    | `str`  | Name/path to the shapefile to plot.                                                         | —       |
| plot_title    | `str`  | Title to assign to plot.                                                                    | `''`    |
| reproject_gcs | `bool` | Determine whether to reproject data before visualization.                                   | `False` |
| projection    | `int`  | If reprojecting data, specifies which EPSG to reproject to.                                 | `4269`  |
| crop_to_shape | `str`  | State shapefile to crop data to.                                                            | `None`  |
| save_to_img   | `str`  | Name/path of file to save plot to.                                                          | `None`  |
| verbose       | `bool` | Determine if additional print statements should be used to track computation of parameters. | `False` |

**Outputs**

| Name  | Type | Description                       | Default |
| ----- | ---- | --------------------------------- | ------- |
| Image | —    | Creates a plot for the shapefile. | —       |

---

### `print_computable_parameters()`

```python
geotiled.print_computable_parameters()
```

Outputs all computable terrain parameters with GEOtiled. A special indicator ‘(G)’ indicates it is computable with GDAL.

---

### `print_region_codes()`

```python
geotiled.print_region_codes()
```

Outputs state abbreviations and their correlating states.

---

### `reproject()`

```python
geotiled.reproject(input_file, output_file, projection, cleanup=False, verbose=False)
```

Reprojects a raster file to a specified projection where the result is saved to a new file. Multithreading is utilized to improve performance.

**Parameters**

| Name        | Type   | Description                                                                                             | Default |
| ----------- | ------ | ------------------------------------------------------------------------------------------------------- | ------- |
| input_file  | `str`  | Name/path of raster to reproject.                                                                       | —       |
| output_file | `str`  | Name/path of file to write reprojected data to.                                                         | —       |
| projection  | `str`  | Projection to use for reprojection. Can be a EPSG code (e.g. EPSG:4269) or the name/path to a WKT file. | —       |
| cleanup     | `bool` | Determine if input_file should be deleted after reprojection is complete.                               | `False` |
| verbose     | `bool` | Determine if additional print statements should be used to track computation.                           | `False` |

> **NOTE:** It is important to use this for preprocessing DEM data to have units of meters.

---

### `set_working_directory()`

```python
geotiled.set_working_directory(path)
```

Creates needed directories then sets the specified path as the working directory. It is important to run this before anything else to ensure data is stored in a desired location.

**Parameters**

| Name | Type  | Description                         | Default |
| ---- | ----- | ----------------------------------- | ------- |
| path | `str` | Working directory to store data in. | —       |

---

### `tif2csv()`

```python
geotiled.tif2csv(input_file, output_file, band_names, verbose=False)
```

Converts a GeoTIFF file and converts it into the CSV format. The CSV columns are the x coordinate, y coordinate, and raster values. NaN values indicate no data found at a particular index. If the GeoTIFF has multiple bands, the user can set the names of subsequent columns in the CSV containing each bands value for each x,y coordinate using the band_names parameter.

**Parameters**

| Name        | Type        | Description                                                                                 | Default |
| ----------- | ----------- | ------------------------------------------------------------------------------------------- | ------- |
| input_file  | `str`       | Name/path to GeoTIFF file to convert to a CSV.                                              | —       |
| output_file | `str`       | Name/path to CSV file to save CSV data to.                                                  | —       |
| band_names  | `List[str]` | Names of columns correlating to raster bands. Order matters.                                | —       |
| verbose     | `bool`      | Determine if additional print statements should be used to track computation of parameters. | `False` |
