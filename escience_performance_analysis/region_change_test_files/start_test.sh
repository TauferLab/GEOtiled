#!/bin/bash

regions=( "flat" "mountain" )
methods=( "GEOtiled-G" "GEOtiled-SG" )
parameters=( "hillshade" "slope" "aspect" "plan_curvature" "profile_curvature" "convergence_index" "filled_depressions" "watershed_basins" "total_catchment_area" "flow_width" "specific_catchment_area" "channel_network" "drainage_basins" "flow_direction" "flow_connectivity" )
runs=( "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" ) 

source_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
configuration_file="$source_dir/configure_test.sh"

for region in ${regions[@]}; do
 for method in ${methods[@]}; do
  for parameter in ${parameters[@]}; do
   for run in ${runs[@]}; do
    $configuration_file "$region" "$method" "$parameter" "$run"
   done
  done
 done
done