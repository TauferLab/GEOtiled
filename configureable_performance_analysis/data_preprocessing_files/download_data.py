import geotiled
import os

try:
    # Download 10m CONUS dataset
    data_10m_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_10m.txt')
    geotiled.download_files(data_10m_file, 'dem_tiles')
    geotiled.build_mosaic('dem_tiles', 'mosaic.tif', cleanup=True)
    geotiled.reproject('mosaic.tif', 'elevation.tif', projection='EPSG:5070', cleanup=True)
except Exception as e:
    print(f"Error downloading files: {e}")

try:
    # Crop test files for different tile size tests
    geotiled.crop_pixels('elevation.tif', 'ts20000.tif', [300000,80000,160000,160000])
    geotiled.crop_pixels('ts20000.tif', 'ts17500.tif', [0,0,140000,140000])
    geotiled.crop_pixels('ts17500.tif', 'ts15000.tif', [0,0,120000,120000])
    geotiled.crop_pixels('ts15000.tif', 'ts12500.tif', [0,0,100000,100000])
    geotiled.crop_pixels('ts12500.tif', 'ts10000.tif', [0,0,80000,80000])
    geotiled.crop_pixels('ts10000.tif', 'ts7500.tif', [0,0,60000,60000])
    geotiled.crop_pixels('ts7500.tif', 'ts5000.tif', [0,0,40000,40000])
    geotiled.crop_pixels('ts5000.tif', 'ts2500.tif', [0,0,20000,20000])
except Exception as e:
    print(f"Error cropping files: {e}")

try:
    # Download and preprocess flat dataset
    geotiled.fetch_dem(bbox={"xmin": -98.00, "ymin": 47.00, "xmax": -96.00, "ymax": 45.00}, dataset="10m", save_to_txt=False, download=True)
    geotiled.build_mosaic('dem_tiles','pre_flat.tif', cleanup=True)
    geotiled.reproject('pre_flat.tif', 'rpre_flat.tif', projection='EPSG:5070', cleanup=True)
    geotiled.crop_pixels('rpre_flat.tif', 'flat.tif', [7000,15000,20000,20000])
except Exception as e:
    print(f"Error generating flat dataset: {e}")

try:
    # Download and preprocess mountain dataset
    geotiled.fetch_dem(bbox={"xmin": -112.00, "ymin": 44.00, "xmax": -110.00, "ymax": 42.00}, dataset="10m", save_to_txt=False, download=True)
    geotiled.build_mosaic('dem_tiles','pre_mountain.tif', cleanup=True)
    geotiled.reproject('pre_mountain.tif', 'rpre_mountain.tif', projection='EPSG:5070', cleanup=True)
    geotiled.crop_pixels('rpre_mountain.tif', 'mountain.tif', [11000,16000,20000,20000])
except Exception as e:
    print(f"Error generating mountain dataset: {e}")