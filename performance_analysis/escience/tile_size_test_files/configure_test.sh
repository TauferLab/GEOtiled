#!/bin/bash

# Input validation
if [ $# -ne 4 ]; then
 echo "Four arguments required."
 exit 1
fi

# Get path to python file to execute
source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python_file="$source_dir/tile_size_test.py"

# Start memory logging
log_name="mem_logs/${2}_${3}_${1}_${4}.csv"
pmrep -t 1sec -b MB -o csv -F $log_name mem.util.used & 
PMVAL_PID=$!

# Begin test
python3 $python_file $@

# End memory logging
kill $PMVAL_PID

# File cleanup
parameter_file="${3}.tif"
rm $parameter_file

if [ "$2" = "GEOtiled-G" ] || [ "$2" = "GEOtiled-SG" ]; then
 parameter_tiles1="${3}_tiles"
 parameter_tiles2="unbuffered_${3}_tiles"
 rm -r $parameter_tiles1
 rm -r $parameter_tiles2
 rm -r elevation_tiles
fi