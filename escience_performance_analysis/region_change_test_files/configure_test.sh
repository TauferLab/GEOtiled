#!/bin/bash

# Input validation
if [ $# -ne 4 ]; then
 echo "Four arguments required."
 exit 1
fi

# Get path to python file to execute
source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python_file="$source_dir/region_change_test.py"

# Start memory logging
log_name="mem_logs/${1}_${2}_${3}_${4}.csv"
pmrep -t 1sec -b MB -o csv -F $log_name mem.util.used & 
PMVAL_PID=$!

# Begin test
python3 $python_file $@

# End memory logging
kill $PMVAL_PID

# File cleanup
if [ "${2}" = "GEOtiled-G" ] && [[ "${3}" = "slope" || "${3}" = "aspect" || "${3}" = "hillshade" ]] || [ "${2}" = "GEOtiled-SG" ]; then
 file_name="${3}.*"
 folder_name1="${3}_tiles"
 folder_name2="unbuffered_${3}_tiles"

 rm $file_name
 rm -r elevation_tiles
 rm -r $folder_name1
 rm -r $folder_name2
fi

if [ "${2}" = "GEOtiled-SG" ] && [ "${3}" = "specific_catchment_area" ]; then
 rm -r total_catchment_area_tiles
fi