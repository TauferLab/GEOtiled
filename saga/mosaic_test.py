import geotiledsaga as gts
import pandas as pd
import time
import os

TILE_SPLIT_SQRT_MIN = 2
TILE_SPLIT_SQRT_MAX = 10

data_path = '/media/volume/geotiled-saga/tile_compute_test_ca'
gts.set_working_directory(data_path)

results = os.path.join(os.getcwd(), 'mosaic_test.csv')
f = open(results, 'w')
f.write('tile_count,execution_time,peak_mem_usage\n')
f.close()

for i in range(TILE_SPLIT_SQRT_MIN,TILE_SPLIT_SQRT_MAX+1):
    # Compute tile count
    tile_count = i**2

    # Set directory to correct tile folder
    gts.set_working_directory(os.path.join(os.getcwd(), str(tile_count) + '_tiles'))

    # Start logging memory usage
    mem_log_file = os.path.join(os.getcwd(), 'mosaic_mem_log.csv')
    f = open(mem_log_file, 'w')
    f.write('mem_usage\n')
    f.close()
    
    start_cmd = ['tmux', 'new-session', '-d', '-s', 'track-mem', 'python', '/home/exouser/GEOtiled/geotiled-saga/track_memory_usage.py', mem_log_file]
    gts.__bash(start_cmd)
    
    # Do mosaicking and record its execution time
    start_time = time.time()
    gts.build_mosaic_buffer('ctt_convergence_index_tiles', 'convergence_index')
    end_time = time.time()

    # End memory logging
    end_cmd = ['tmux','kill-session','-t','track-mem']
    gts.__bash(end_cmd)

    # Compute peak memory usage
    mem_data = pd.read_csv(mem_log_file)
    df = pd.DataFrame(mem_data)
    mem_usage = df['mem_usage']
    peak_usage = max(mem_usage) - min(mem_usage)

    # Write results to file
    formatted_results = str(tile_count) + ',' + str(end_time-start_time) + ',' + str(peak_usage) + '\n'
    f = open(results, 'a')
    f.write(formatted_results)
    f.close()

    # Reset working directory
    gts.set_working_directory(data_path)