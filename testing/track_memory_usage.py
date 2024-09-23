import psutil
import time
import sys

def log_mem_usage(log_file, record_time=0.1):
    while(1):
        current_memory_usage = psutil.virtual_memory().total - psutil.virtual_memory().available
        formatted_result = str(current_memory_usage) + '\n'
        f = open(log_file, 'a')
        f.write(formatted_result)
        f.close()
        time.sleep(record_time)

log_file_path = sys.argv[1]
log_mem_usage(log_file_path)