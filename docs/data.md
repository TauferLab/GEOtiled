# 💾 Data Guide

The GEOtiled workflow is designed to work with files in the [GeoTIFF](https://www.earthdata.nasa.gov/about/esdis/esco/standards-practices/geotiff) (extension: `.tif`) format, and allows for download of [Digital Elevation Model (DEM)](https://www.usgs.gov/faqs/what-a-digital-elevation-model-dem) data from the [United States Geological Survey (USGS)](https://www.usgs.gov/). This page goes over available data from the USGS useable by GEOtiled, the list of computable terrain parameters, and processing recommendations to ensure computational compatibility and quality.

If [GDAL](https://gdal.org/en/stable/) is installed, it is possible to view metadata of a GeoTIFF file such as the total size (x,y), coordinate reference system (CRS), pixel size, coordinate extents, and band(s) information using the following command line function `gdalinfo <your_file>.tif`.

## Compatible Datasets

The USGS offers DEMs at various resolutions on the [TNM Download](https://apps.nationalmap.gov/downloader/#/) page. GEOtiled always downloads the **latest versions** of available DEMs from the USGS. Below are the available datasets compatible with GEOtiled and along with their resolution and spatial coverage.

National Elevation Dataset (NED) Alaska 2 arc-second Current  
&nbsp;&nbsp; \- Resolution: 60 meters  
&nbsp;&nbsp; \- Coverage: Alaska  
National Elevation Dataset (NED) 1 arc-second Current  
&nbsp;&nbsp; \- Resolution: 30 meters  
&nbsp;&nbsp; \- Coverage: USA, Canada, Mexico  
National Elevation Dataset (NED) 1/3 arc-second Current  
&nbsp;&nbsp; \- Resolution: 10 meters  
&nbsp;&nbsp; \- Coverage: Continental United States (CONUS)  
Alaska IFSAR 5 meter DEM  
&nbsp;&nbsp; \- Resolution : 5 meters  
&nbsp;&nbsp; \- Coverage: Alaska  
Digital Elevation Model (DEM) 1 meter  
&nbsp;&nbsp; \- Resolution: 1 meter  
&nbsp;&nbsp; \- Coverage: CONUS (partial)

## Computable Terrain Parameters

The following is a comprehensive list of computable terrain parameters with GEOtiled using either underlying GIS library. SAGA documentation for each terrain parameter's underlying function call is linked. All terrain parameters computed by GDAL are done using the [gdal.DEMProcessing()](https://gdal.org/en/stable/api/python/utilities.html) function.

### GDAL and SAGA

[Slope](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_0.html)  
[Aspect](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_0.html)  
[Hillshade](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_lighting_0.html)  
&nbsp;&nbsp; \- GDAL and SAGA define hillshade differently, so results will vary substantially between the two

### SAGA Only

[Plan Curvature](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_0.html)  
[Profile Curvature](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_0.html)  
[Convergence Index](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_1.html)  
[Total Catchment Area](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_hydrology_0.html)  
[Flow Width](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_hydrology_19.html)  
[Specific Catchment Area](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_hydrology_19.html)  
&nbsp;&nbsp; \- Total Catchment Area is a required input for computation, and is computed automatically by GEOtiled  
[Channel Network (Shapefile)](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_5.html)  
[Channel Network (GeoTIFF)](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_5.html)  
[Drainage Basins (Shapefile)](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_5.html)  
[Drainge Basins (GeoTIFF)](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_5.html)  
[Flow Connectivity](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_5.html)  
[Flow Direction](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_5.html)  
[Channel Network Base Level](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_3.html)  
&nbsp;&nbsp; \- Channel Network (GeoTIFF) is a required input for computation, and is computed automatically by GEOtiled  
[Channel Network Distance](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_channels_3.html)  
&nbsp;&nbsp; \- Channel Network (GeoTIFF) is a required input for computation, and is computed automatically by GEOtiled  
[Filled Depressions](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_preprocessor_4.html)  
&nbsp;&nbsp; \- It is recommended to compute this with a DEM that has a coordinate reference system (CRS) in units of degrees with square grid cells  
[Watershed Basins](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_preprocessor_4.html)  
&nbsp;&nbsp; \- It is recommended to compute this with a DEM that has a coordinate reference system (CRS) in units of degrees with square grid cells  
[LS Factor](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_hydrology_22.html)  
&nbsp;&nbsp; \- Both Slope and Specific Catchment Area are required inputs for computation, and they are automatically computed by GEOtiled  
[Topographic Wetness Index](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_hydrology_20.html)  
&nbsp;&nbsp; \- Both Slope and Specific Catchment Area are required inputs for computation, and they are automatically computed by GEOtiled  
[Valley Depth](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_14.html)  
&nbsp;&nbsp; \- The underlying function for this terrain parameter runs much slower than other functions.  
[Relative Slope Position](https://saga-gis.sourceforge.io/saga_tool_doc/9.3.1/ta_morphometry_14.html)  
&nbsp;&nbsp; \- The underlying function for this terrain parameter runs much slower than other functions.  

## Processing Recommendations

### Coordinate Reference System (CRS)

DEMs from the USGS in the GeoTIFF format are stored in the [EPSG:4269](https://epsg.io/4269) coordinate system. The CRS fully covers North America and has dimensions of degrees. It is recommended that the projection to use when computing terrain parameters with GEOtiled:

* Has pixel sizes with even length and width (i.e., square) as it is a computational requirement for underlying GIS libraries GDAL and SAGA 
* Has full coverage of at least the Continental United States (CONUS) to ensure the accuracy of the final data
* Has units of meters to prevent specific computational issues for most terrain parameters

A recommended projection to use is the [EPSG:5070](https://epsg.io/5070) CRS. It covers North America, and has units of meters with square grid cells.

GeoTIFF files can be reprojected using the `geotiled.reproject()` function (whose details are in the [Library Guide](./library.md)).

**IMPORTANT**: Reprojection will cause changes to the size (x,y) of the data. Ensure that other configurations set within GEOtiled (such as tile size) are set accordingly to account for changes in input data caused by reprojection.

### Setting Tile Size

For the best performance results during computation of terrain parameters, the number of tiles an input DEM is split into should match the number of cores being used for multiprocessing. You can use `gdalinfo <your_file>.tif` in the command line to look at the size (x,y) of the input file to determine the best tile split to use. Optionally, the `num_tiles` parameter in the `geotiled.crop_and_compute()` function (whose details are in the [Library Guide](./library.md)) can be set to the number of cores to automatically determine good tile size dimensions. 
