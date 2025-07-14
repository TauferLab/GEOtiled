# Script implementing functions to get accuracy metrics of two-dimensional raster data 
# Last Updated: 12/23/2024
# Author(s): Gabriel Laboy (@glaboy-vol)

from osgeo import gdal
import numpy as np

def absolute_difference(file1, file2):
    """
    Calculates absolute difference.

    Computes absolute difference of two raster files that have
    data in the form of two-dimensional arrays.

    Parameters
    ----------
    file1 : str
        Path to first file to compare.
    file2 : str
        Path to second file to compare.
    """

    # Open and convert raster data to two-dimensional arrays
    data1 = gdal.Open(file1)
    array1 = data1.ReadAsArray()
    data2 = gdal.Open(file2)
    array2 = data2.ReadAsArray()
    
    # Cast two-dimensional arrays to numpy arrays
    np_array1, np_array2 = map(np.array, (array1, array2))

    # Sum together absolute differences of each array value to get absolute difference
    absolute_difference = np.sum(np.absolute(np_array1 - np_array2))
    return absolute_difference

def mse(file1, file2):
    """
    Calculates mean squared error.

    Computes mean squared error of two raster files that have
    data in the form of two-dimensional arrays.

    Parameters
    ----------
    file1 : str
        Path to first file to compare.
    file2 : str
        Path to second file to compare.
    """

    # Open and convert raster data to two-dimensional arrays
    data1 = gdal.Open(file1)
    array1 = data1.ReadAsArray()
    data2 = gdal.Open(file2)
    array2 = data2.ReadAsArray()

    # Cast two-dimensional arrays to numpy arrays
    np_array1, np_array2 = map(np.array, (array1, array2))

    # Average of differences of array values squared to get mse
    mse = np.mean((np_array1 - np_array2) ** 2)
    return mse