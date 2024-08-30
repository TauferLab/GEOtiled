from osgeo import gdal
import numpy as np

import psutil
import sys
import os

gdal.UseExceptions()

file = '/media/volume/geotiled-saga/tile_compute_test_de/elevation.tif'

file_size = os.path.getsize(file)
matrix_data = gdal.Open(file)
array_data = matrix_data.ReadAsArray()
nodata_percentage = (len(array_data[array_data == -999999.0]) / (len(array_data) * len(array_data[0]))) * 100
print("%s | FS: %i B | NDP: %.2f" % ([len(array_data[0]), len(array_data)], file_size, nodata_percentage))

# def printme():
#     print('coolstuff')

# with open('file.txt', 'a') as sys.stdout:
#     print('Current RAM Usage: %.3f MB' % ((psutil.virtual_memory().total - psutil.virtual_memory().available) / (1024*1024)))
#     printme()

# print('Current RAM Usage: %.3f MB' % ((psutil.virtual_memory().total - psutil.virtual_memory().available) / (1024*1024)))