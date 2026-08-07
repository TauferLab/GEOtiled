# 📦 GEOtiled Library

The `geotiled` library offers numerous functions for download and preprocessing of Digital Elevation Model (DEM) data and computation of terrain parameters from said DEMs. Once you have cloned the the GitHub repository, this will guide you through available functions from the library along with a quick start guide on core functions to use for generating terrain parameters.

## 🚀 Quick Start

Before installation, ensure you have [Python](https://www.python.org/downloads/) 3.7.0 or higher and [pip](https://pip.pypa.io/en/stable/installation/) installed. GEOtiled was developed and tested in Ubuntu 24.04, and stability on other Linux distributions remains untested.

### Environment Setup

1. Create a new Python virtual environment in a desired directory
    > Note: `<your_path>` should be replaced with your desired working directory
```
python3 -m venv <your_path>/geotiled_env
```
2. Activate the new virtual environment
```
source <your_path>/geotiled_env/bin/activate
```
3. Update pip
```
pip install --upgrade pip
```

### System Packages

1. Update apt-get
```
sudo apt-get update
```
2. Install libgdal-dev
```
sudo apt install libgdal-dev=3.8.4+dfsg-3ubuntu3
```
3. Install SAGA
```
sudo apt-get install saga
```
4. [**Optional**] Install Performance Copilot (Only used for performance analysis notebooks)
```
sudo apt-get install pcp-zeroconf
```

### Installing GEOtiled

1. Ensure virutal environment is activated
```
source <your_path>/geotiled_env/bin/activate
```
2. Clone the repository in a desired working directory
```
git clone https://github.com/TauferLab/GEOtiled
```
3. Change to the geotiled directory
```
cd <your_path>/GEOtiled/geotiled
```
4. Install editable library
```
pip install -e .
```

### Common Issues

If the following error occurs during computation:

`ImportError: cannot import name '_gdal_array' from 'osgeo'`

Then pip installed a cached version of GDAL. Run the following in the virtual environment to correct the issue:
```
pip install --no-cache --force-reinstall gdal[numpy]==3.8.4
```

## 📚 GEOtiled Library

### Getting Started

The `geotiled` library can be imported by adding the following line to the top of a Python script:

```python
import geotiled
```

It is recommended to set a working directory for storing data, as GEOtiled produces a lot of intermediary and final files during computation. GEOtiled comes with a function to set the working directory (which also handles creating missing folders).

```python
geotiled.set_working_directory("your/working/directory")
```

### Downloading DEMs

`geotiled.fetch_dems()` downloads DEMs directly from the USGS [TNM Download](https://apps.nationalmap.gov/downloader/#/) page within the bounds of a specified state or coordinate extents along with a desired resolution. Available resolutions for download can be viewed under the *Compatible Datasets* section in the [Data Guide](./data.md). The user has the option to either save the USUS download URLs to a text file and/or download the DEMs immediately. 

`geotiled.download_files()` downloads files from a list of URLs provided in a text file. The downloaded files will be stored in a user-specified output folder located in the working directory.

An example to download 30 meter DEMs for the state of Tennessee is given below:

```python
geotiled.fetch_dems(shapefile="TN", dataset="30m", txt_file="TN_30m_urls.txt", download=False)

geotiled.download_files(download_list="TN_30m_urls.txt", download_folder="dem_tiles")
```

`geotiled.mosaic_rasters()` merges together numerous raster files into a single raster file of one band. It is important to merge together DEMs from the USGS before computation of terrain parameters. It will merge together all GeoTIFF files in a specified folder.

`geotiled.reproject()` changes the coordinate reference system (CRS) of raster data. DEMs from the USGS come in the [EPSG:4269](https://epsg.io/4269) CRS which uses units of decimal degrees. In general, it is good practice to reproject input data into a CRS with units of meters and square grid cells, so DEMs from the USGS should be reprojected before computation. The function allows for reprojection using a valid EPSG code or a [WKT](https://libgeos.org/specifications/wkt/) file with projection information.

An example for preprocessing DEM data is given below (using CRS [EPSG:5070](https://epsg.io/5070)):

```python
geotiled.mosaic_rasters(input_folder="dem_tiles", output_file="mosaic.tif")

geotiled.reproject(input_file="mosaic.tif", output_file="elevation.tif", projection="EPSG:5070")
```

### Computing Terrain Parameters

`geotiled.crop_and_compute()` handles the cropping of the input DEM and concurrent computation of terrain parameters from subsequent DEM tiles. The user must specify the list of terrain parameters to compute. The full list of computable terrain parameters can be found under the *Computable Terrain Parameters* section in the [Data Guide](./data.md). Additionally, the user must also specify the tile size for each cropped tile and which library to use for computation. The function automatically consumes the max number of CPU cores on a machine for concurrent computation unless a specific count is specified.

`geotiled.mosaic_rasters()` is again used for merging terrain parameter tiles in the GeoTIFF format. `geotiled.crop_and_compute()` also handles debuffering of computed terrain parameter tiles, so the correct folder where the terrain parameter tiles are located must be given for accurate final results.

`geotiled.merge_shapefiles()` is used for merginging terrain parameter tiles that are output in the shapefile format. As with `geotiled.mosaic_rasters()`, the user must specify a folder where all shapefiles to merge are located.

An example for computing terrain parameters is given below:

```python
geotiled.crop_and_compute(input_file="elevation.tif", parameter_list=["slope","channel_network"], tile_dimensions=[9103,4195], compute_method="SAGA")

geotiled.mosaic_rasters(input_folder="unbuffered_slope_tiles", output_file="slope.tif")

geotiled.merge_shapefiles(input_folder="channel_network_tiles", output_file="channel_network.shp")
```

### North American Ecoregion Workflow

EPA North American ecoregion codes can be used as DEM sources. GEOtiled infers Level I, II, or III from codes such as 11, 11.1, and 11.1.1, respectively. The EPA boundary archive for the inferred level is downloaded once and cached under the working directory's shapefiles folder.

The standard ecoregion bundle contains 16 terrain parameters: 14 raster products plus channel-network and drainage-basin vectors. Together with the cropped elevation raster, the workflow creates 17 deliverables.

Compute terrain parameters from the complete DEM mosaic before applying the ecoregion cutline. This retains neighboring elevation cells at the boundary while SAGA computes terrain derivatives.

~~~python
region = "11.1.1"

geotiled.fetch_dems(
    ecoregion=region,
    dataset="30m",
    txt_file="ecoregion_30m_urls.txt",
)
geotiled.download_files("ecoregion_30m_urls.txt", "dem_tiles")
geotiled.mosaic_rasters("dem_tiles", "dem_mosaic.tif")
geotiled.reproject("dem_mosaic.tif", "elevation_source.tif", "EPSG:5070")

geotiled.crop_and_compute(
    input_file="elevation_source.tif",
    parameter_list=geotiled.ECOREGION_TERRAIN_PARAMETERS,
    num_tiles=8,
    compute_method="SAGA",
    projection=5070,
)

geotiled.crop_to_ecoregion("elevation_source.tif", "elevation.tif", region)

vector_parameters = {"channel_network", "drainage_basins"}
for parameter in geotiled.ECOREGION_TERRAIN_PARAMETERS:
    if parameter in vector_parameters:
        geotiled.merge_shapefiles(
            input_folder=f"{parameter}_tiles",
            output_file=f"{parameter}.zip",
            ecoregion=region,
        )
    else:
        source = f"{parameter}_source.tif"
        geotiled.mosaic_rasters(f"unbuffered_{parameter}_tiles", source)
        geotiled.crop_to_ecoregion(source, f"{parameter}.tif", region)
~~~

Level I regions can intersect many DEM tiles. Review the generated URL list before enabling immediate downloads.

To generate Level I regions `1` through `15` as separate sequential jobs, use the [Level I command-line workflow](./functions.md#level-i-command-line-workflow). The documented driver accepts one region code at a time:

```bash
python run_level1_ecoregion.py 1 --output-root ./geotiled_level1
python run_level1_ecoregion.py 2 --output-root ./geotiled_level1
```

All 15 can be processed in order with a Bash loop:

```bash
for region in {1..15}; do
    python run_level1_ecoregion.py "$region" --output-root ./geotiled_level1 || break
done
```


### Visualizing Data

`geotiled.plot_raster()` visualizes a GeoTIFF. All plotting features can be found in the [Function Documentation](./functions.md)

`geotiled.plot_shapefile()` visualizes a shapefile. All plotting features can be found in the [Function Documentation](./functions.md)

Examples of how to plot both types of data are given below:

```python
geotiled.plot_raster(input_file="slope.tif", plot_title="TN Slope 30m", reproject_gcs=True, remove_nans=True, shapefiles=["TN"], downsample=5, zunit="Radian", xyunit="Degree", ztype="Slope")

geotiled.plot_shapefile(input_file="channel_network.shp", plot_title="TN Channel Network 30m", reproject_gcs=True, crop_to_shape="TN")
```