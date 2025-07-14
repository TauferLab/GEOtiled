#!/bin/bash

# Input validation
if [ $# -ne 7 ]; then
 echo "Seven arguments required."
 exit 1
fi

# Get path to python file to execute
source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python_file="$source_dir/optimization_test.py"

# Start memory logging
log_name="mem_logs/${1}_${4}_${5}_${3}_${7}.csv"
pmrep -t 1sec -b MB -o csv -F $log_name mem.util.used & 
PMVAL_PID=$!

# Begin test
python3 $python_file $@

# End memory logging
kill $PMVAL_PID

# File cleanup
output_file="${5}.tif"
output_tiles="${5}_tiles"

rm $output_file
rm -r $output_tiles
rm -r elevation_tiles

if [ "${1}" = "optimized" ]; then
 rm -r unbuffered_tiles
fi