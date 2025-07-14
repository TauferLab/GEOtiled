#!/bin/bash

methods=( "GEOtiled-G" "GEOtiled-SG" )
processes=( "1" "2" "4" "8" "16" "32" "64" )
parameters=( "slope" "aspect" "hillshade" )
runs=( "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" ) 

source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
configuration_file="$source_dir/configure_test.sh"

for method in ${methods[@]}; do
 for process in ${processes[@]}; do
  for parameter in ${parameters[@]}; do
   for run in ${runs[@]}; do
    $configuration_file "$method" "$process" "$parameter" "$run"
   done
  done
 done
done