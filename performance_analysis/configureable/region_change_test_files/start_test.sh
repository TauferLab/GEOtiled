#!/bin/bash

# Input validation
if [ $# -ne 6 ]; then
 echo "Six arguments required."
 exit 1
fi

# Get path to python file to execute
source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python_file="$source_dir/region_change_test.py"

# Start memory logging
log_name="mem_logs/${2}_${1}_${4}_${3}_${6}.csv"
pmrep -t 1sec -b MB -o csv -F $log_name mem.util.used & 
PMVAL_PID=$!

# Begin test
python3 $python_file $@

# End memory logging
kill $PMVAL_PID

# File cleanup
if [ "${1}" = "GEOtiled-G" ] && [[ "${4}" = "slope" || "${4}" = "aspect" || "${4}" = "hillshade" ]] || [ "${1}" = "GEOtiled-SG" ]; then
 file_name="${4}.*"
 folder_name1="${4}_tiles"
 folder_name2="unbuffered_${4}_tiles"

 rm $file_name
 rm -r elevation_tiles
 rm -r $folder_name1
 rm -r $folder_name2
fi

if [ "${1}" = "GEOtiled-SG" ] && [ "${4}" = "specific_catchment_area" ]; then
 rm -r total_catchment_area_tiles
fi