#!/bin/bash

# Input validation
if [ $# -ne 3 ]; then
 echo "Three arguments required."
 exit 1
fi

# Get path to python file to execute
source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python_file="$source_dir/optimization_test.py"

# Start memory logging
log_name="mem_logs/${2}_${1}_${3}.csv"
pmrep -t 1sec -b MB -o csv -F $log_name mem.util.used & 
PMVAL_PID=$!

# Begin test
python3 $python_file $@

# End memory logging
kill $PMVAL_PID

# File cleanup
rm slope.tif
rm -r slope_tiles
rm -r elevation_tiles

if [ "${2}" = "optimized" ]; then
 rm -r unbuffered_tiles
fi