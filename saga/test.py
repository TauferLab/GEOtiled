from osgeo import osr, ogr, gdal
from pathlib import Path
from tqdm import tqdm

import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import numpy as np

import concurrent.futures
import multiprocessing
import subprocess
import requests
import tempfile
import zipfile
import shutil
import psutil
import math
import glob
import time
import json
import os
import re

COMPUTABLE_PARAMETERS = ["slope", "aspect", "hillshade", "plan_curvature", "profile_curvature", "convergence_index"]

def determine_if_path(path):
    """
    Determine if entered path is a folder or filename or full path.

    Determines if an entered path is a folder or filename or full path.
    If the entered path is a full, valid path, the path is returned.
    If the entered path is only a folder or filename, a path is built 
    with the specified name and current working directory and the 
    new full path is returned.

    Parameters
    ----------
    path : str
        A full path or folder/file name.
    """

    # Determine if the entered path is just a folder/filename 
    if Path(path).exists():
        return path
    else:
        new_path = os.path.join(os.getcwd(), path)
        Path(new_path).mkdir(parents=True, exist_ok=True)
        return new_path

print(determine_if_path('/media/geotiled/geotiled-saga/tile_compute_test_tn'))

print(determine_if_path('test_folder'))

print(determine_if_path('test.txt'))