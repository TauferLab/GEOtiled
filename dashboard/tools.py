###############
### IMPORTS ###
###############

import os
import yaml
import math
import logging
import zipfile
import warnings

from botocore.client import Config
from boto3.session import Session
from dotenv import load_dotenv

import matplotlib.pyplot as plt
from osgeo import gdal, osr
import geopandas as gpd
import numpy as np

logger = logging.getLogger(__name__)
logging.basicConfig(filename="debug.log", encoding="utf-8", level=logging.INFO)

# Enable exceptions for GDAL 
gdal.UseExceptions()

# Suppress specific warning related to loading shapefiles
warnings.filterwarnings("ignore", message=r"Measured \(M\) geometry types are not supported.*")

#################
### FUNCTIONS ###
#################

def get_aws_bucket(verify=True):
    """
    Load AWS bucket given config
    """
    load_dotenv()

    config = Config(signature_version="s3v4")
    endpoint_url = os.getenv("ENDPOINT_URL")
    bucket_name = os.getenv("BUCKET_NAME")
    aws_access_key_id = os.getenv("AWS_ACCESS_KEY_ID")
    aws_secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")

    bucket = (
        Session()
        .resource(
            "s3",
            endpoint_url=endpoint_url,
            config=config,
            verify=verify,
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
        )
        .Bucket(bucket_name)
    )
    return bucket


def download(resolution, state, terrain_parameter, storage_directory='./'):
    # Create the directory to store files
    full_storage_directory = os.path.join(storage_directory,f"utk/conus/{resolution}/{state}")
    os.makedirs(full_storage_directory, exist_ok=True)

    # Load the configuration file
    config_path = os.getenv("CONFIG_PATH", "./config.yaml")
    config = {}
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"could not initialize configuration, configuration path does not exists: {e}")
        return

    # Download the requested terrain parameter
    file = os.path.join(
        config["prefix"], config["resolutions"][resolution], config["states"][state], config["files"][terrain_parameter]
    )

    # destination (your host path)
    dst = os.path.join(storage_directory, file)
    dst_dir = os.path.dirname(dst)
    os.makedirs(dst_dir, exist_ok=True)

    if os.path.exists(dst):
        print(f"{dst} already exists. Skipping download.")
        return
    
    try:
        s3 = get_aws_bucket()
        s3.download_file(file, dst)
        logger.info(f"downloaded {file}")
    except Exception as e:
        logger.error(e)

    try:
        # Unzip any shapefiles
        if (terrain_parameter == 'CHANNEL_NETWORK') or (terrain_parameter == 'DRAINAGE_BASINS'):
            with zipfile.ZipFile(dst, 'r') as file:
                file.extract(f"{state}/{terrain_parameter.lower()}.dbf", dst_dir)
                file.extract(f"{state}/{terrain_parameter.lower()}.shp", dst_dir)
                file.extract(f"{state}/{terrain_parameter.lower()}.shx", dst_dir)

            os.rename(os.path.join(dst_dir,f"{state}/{terrain_parameter.lower()}.dbf"), os.path.join(dst_dir,f"{terrain_parameter.lower()}.dbf"))
            os.rename(os.path.join(dst_dir,f"{state}/{terrain_parameter.lower()}.shp"), os.path.join(dst_dir,f"{terrain_parameter.lower()}.shp"))
            os.rename(os.path.join(dst_dir,f"{state}/{terrain_parameter.lower()}.shx"), os.path.join(dst_dir,f"{terrain_parameter.lower()}.shx"))
            os.rmdir(os.path.join(dst_dir,state))
    except Exception as e:
        logger.error(e)


def visualize(resolution, state, terrain_parameter, downsample=1, cmap="inferno", dpi=300, vmin=None, vmax=None, nancolor="white"):
    # Load the configuration file
    config_path = os.getenv("CONFIG_PATH", f"./config.yaml")
    config = {}
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(
            f"could not initialize configuration, configuration path does not exists: {e}"
        )
        return
    
    # The file to load
    file = os.path.join(
        config["prefix"], config["resolutions"][resolution], config["states"][state], config["files"][terrain_parameter]
    )
    file_path = os.path.join('./', file)

    # Visualize data
    try:
        if ".tif" in file: 
            # Get raster metadata
            dataset = gdal.Open(file_path)
            band = dataset.GetRasterBand(1)
            geotransform = dataset.GetGeoTransform()
            spatial_ref = osr.SpatialReference(wkt=dataset.GetProjection())
        
            # Extract spatial information about raster
            proj_name = spatial_ref.GetAttrValue("PROJECTION")
            proj_name = proj_name if proj_name else "GCS, No Projection"
            data_unit = "" or spatial_ref.GetLinearUnitsName()
            coord_unit = "" or spatial_ref.GetAngularUnitsName()
            z_type = "" if band.GetDescription() == "" else band.GetDescription()
    
            # Downsampled data
            raster_array = gdal.Warp("", file_path, format="MEM", 
                                     width=int(dataset.RasterXSize/downsample), 
                                     height=int(dataset.RasterYSize/downsample)).ReadAsArray()
        
            # Mask nodata values
            raster_array = np.ma.array(raster_array, mask=np.equal(raster_array, band.GetNoDataValue()))
    
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
    
            # Adjust colorbar and title
            cbar = fig.colorbar(sm, fraction=0.046*raster_array.shape[0]/raster_array.shape[1], pad=0.04)
            cbar_ticks = np.linspace(vmn, vmx, 8)
            cbar.set_ticks(cbar_ticks)
            if data_unit == "unknown":
                if terrain_parameter in ['ELEVATION','CLOSED_DEPRESSIONS','DRAINAGE_BASINS_GRID','FLOW_WIDTH']:
                    cbar.set_label("Meters")
                elif terrain_parameter in ['SLOPE','ASPECT']:
                    cbar.set_label("Radians")
                elif terrain_parameter == 'HILLSHADE':
                    cbar.set_label("Shade Level")
                elif terrain_parameter in ['PLAN_CURVATURE','PROFILE_CURVATURE']:
                    cbar.set_label(r"Meters$^{-1}$")
                elif terrain_parameter == 'TOTAL_CATCHMENT_AREA':
                    cbar.set_label(r"Meters$^{2}$")
                elif terrain_parameter in ['CHANNEL_NETWORK_GRID', 'FLOW_CONNECTIVITY']:
                    cbar.set_label("Presence")
                elif terrain_parameter == 'FLOW_DIRECTION':
                    cbar.set_label("Cardinal Direction (0 is north)")
                elif terrain_parameter == 'CONVERGENCE_INDEX':
                    cbar.set_label("Convergence Index")
                elif terrain_parameter == 'WATERSHED_BASINS':
                    cbar.set_label(r"Kilometers$^{2}$")
                else:
                    cbar.set_label("")
            else:
                cbar.set_label(f"{z_type} ({data_unit})")
    
            ax.set_title(f"{terrain_parameter.replace('_', ' ')} for {state} {resolution}", fontsize=12, fontweight="bold", pad=20)
            ax.tick_params(axis="both", which="both", bottom=True, top=False, left=True, right=False, color="black", length=5, width=1)
            
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
            y_label = f"Latitude ({coord_unit})" if spatial_ref.EPSGTreatsAsLatLong() else f"Northing ({coord_unit})"
            x_label = f"Longitude ({coord_unit})" if spatial_ref.EPSGTreatsAsLatLong() else f"Easting ({coord_unit})"
            ax.set_ylabel(y_label)
            ax.set_xlabel(x_label)
        
            # ax.ticklabel_format(style="plain", axis="both")  # Prevent scientific notation on tick labels
            ax.set_aspect("equal")

            return fig
        elif ".zip" in file:
            file = file.replace('.zip','.shp')

            # Load the shapefile
            input_data = gpd.read_file(file)
            input_data.set_crs(epsg=4269)
            
            # Get bounds of shapefile
            minx, miny, maxx, maxy = input_data.total_bounds
            
            # Set up the ticks for x and y axis
            x_ticks = np.linspace(minx, maxx, 5)
            y_ticks = np.linspace(miny, maxy, 5)
        
            # Format the tick labels to two decimal places
            x_tick_labels = [f"{tick:.2f}" for tick in x_ticks]
            y_tick_labels = [f"{tick:.2f}" for tick in y_ticks]
            
            # Plot the shapefile
            fig, ax = plt.subplots()
            input_data.plot(ax=ax)
                
            plt.title(f"{terrain_parameter.replace('_', ' ')} for {state} {resolution}", fontsize=12, fontweight='bold', pad=20)
            plt.xlabel("Longitude (Degrees)", fontsize=10)
            plt.ylabel("Latitude (Degrees)", fontsize=10)
            ax.set_xlim(minx, maxx)
            ax.set_ylim(miny, maxy)
            ax.set_xticks(x_ticks)
            ax.set_yticks(y_ticks)
            ax.set_xticklabels(x_tick_labels, fontsize=9)
            ax.set_yticklabels(y_tick_labels, fontsize=9)
            
            ax.set_aspect("equal")

            return fig
    except Exception as e:
        logger.error(e)