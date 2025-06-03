#!/bin/bash

source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python_file="$source_dir/download_data.py"

python3 $python_file