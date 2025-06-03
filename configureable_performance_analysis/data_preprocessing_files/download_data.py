from osgeo import gdal
import geotiled
import os

def crop(input_file, output_file, window, cleanup=False):
    # Window to crop by [left_x, top_y, width, height]
    translate_options = gdal.TranslateOptions(srcWin=window, creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"])
    
    # Perform translation
    gdal.Translate(output_file, input_file, options=translate_options)

    # Remove input file if requested
    if cleanup:
        os.remove(input_file)

# Download 10m CONUS dataset
data_10m_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_10m.txt')
geotiled.download_files(data_10m_file, 'dem_tiles')
geotiled.build_mosaic('dem_tiles','mosaic.tif', cleanup=True)

# Crop initial file for different tile size tests
crop('mosaic.tif', 'ts20000.tif', [90000,30000,160000,160000])
crop('ts20000.tif', 'ts17500.tif', [0,0,140000,140000])
crop('ts20000.tif', 'ts15000.tif', [0,0,120000,120000])
crop('ts20000.tif', 'ts12500.tif', [0,0,100000,100000])
crop('ts20000.tif', 'ts10000.tif', [0,0,80000,80000])
crop('ts20000.tif', 'ts7500.tif', [0,0,60000,60000])
crop('ts20000.tif', 'ts5000.tif', [0,0,40000,40000])
crop('ts20000.tif', 'ts2500.tif', [0,0,20000,20000])

# Download and preprocess flat dataset
geotiled.fetch_dem(bbox={"xmin": -98.00, "ymin": 47.00, "xmax": -96.00, "ymax": 45.00}, dataset="10m", save_to_txt=False, download=True)
geotiled.build_mosaic('dem_tiles','pre_flat.tif', cleanup=True)
crop('pre_flat.tif', 'flat.tif', [0, 0, 20000, 20000], cleanup=True)

# Download and preprocess mountain dataset
geotiled.fetch_dem(bbox={"xmin": -112.00, "ymin": 44.00, "xmax": -110.00, "ymax": 42.00}, dataset="10m", save_to_txt=False, download=True)
geotiled.build_mosaic('dem_tiles','pre_mountain.tif', cleanup=True)
crop('pre_mountain.tif', 'mountain.tif', [0, 0, 20000, 20000], cleanup=True)