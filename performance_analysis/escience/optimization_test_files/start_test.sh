#!/bin/bash

tile_sizes=( "2500" "15000" )
methods=( "unoptimized" "optimized" )
runs=( "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" ) 

source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
configuration_file="$source_dir/configure_test.sh"

for tile_size in ${tile_sizes[@]}; do
 for method in ${methods[@]}; do
  for run in ${runs[@]}; do
   $configuration_file "$tile_size" "$method" "$run"
  done
 done
done