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

def write_results_to_csv(results_file, results):
    formatted_results = ''
    for result in results:
        formatted_results += str(result) + ','
    formatted_results = formatted_results[:-1] + '\n'
    f = open(results_file, 'a')
    f.write(formatted_results)
    f.close()

# Compute and save averages for the single compute test to a seperate file
averages_file = '/media/volume/geotiled-saga/tile_compute_test_full/single_compute_averages_test.csv'
f = open(averages_file, 'w')
f.write('tile_count,parameter,avg_ex_time,std_ex_time,avg_peak_mem_usage,std_peak_mem_usage\n')
f.close()

for tile_count in [4,9]:
    data = pd.read_csv('/media/volume/geotiled-saga/tile_compute_test_full/single_compute_test.csv')
    df = pd.DataFrame(data)
    df = df[df['tile_count'] == tile_count]
    
    for parameter in ['hillshade','slope']:
        ex_times = df[df['parameter'] == parameter]['ex_time']
        mem_usages = df[df['parameter'] == parameter]['peak_mem_usage']

        ex_time_avg = ex_times.mean()
        ex_time_std = ex_times.std()
        mem_usage_avg = mem_usages.mean()
        mem_usage_std = mem_usages.std()
    
        write_results_to_csv(averages_file, [tile_count, parameter, ex_time_avg, ex_time_std, mem_usage_avg, mem_usage_std])