###############
### IMPORTS ###
###############

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes, mark_inset
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from tabulate import tabulate
from pathlib import Path
from osgeo import gdal
import pandas as pd
import numpy as np
import subprocess
import geotiled
import ast
import os

gdal.UseExceptions()

PARAMETER_CODES = {"hillshade":"H", "slope": "S", "aspect": "A", "plan_curvature": "PlC", "profile_curvature": "PrC", 
                   "convergence_index": "CI", "filled_depressions": "FiD", "watershed_basins": "WB", "total_catchment_area": "TCA", 
                   "flow_width": "FW", "specific_catchment_area": "SCA", "channel_network": "CN", "drainage_basins": "DB", 
                   "flow_direction": "FlD", "flow_connectivity": "FC"}

PARAMETER_GROUPINGS = {"hillshade":["hillshade"], 
                       "curvature":["slope","aspect","plan_curvature","plan_curvature"],
                       "convergence_index":["convergence_index"],
                       "depressions":["filled_depressions","watershed_basins"],
                       "total_catchment_area":["total_catchment_area"],
                       "specific_catchment_area":["flow_width","specific_catchment_area"],
                       "channel_network":["channel_network","drainage_basins","flow_direction","flow_connectivity"]}

############################
### PROCESSING FUNCTIONS ###
############################

def get_file_directory():
    return os.path.dirname(os.path.abspath(__file__))


def get_peak_memory_usage(mem_csv, start_index, end_index):
    df_mem = pd.read_csv(mem_csv)
    df_filt = df_mem[start_index:end_index+1]
    mu_max = df_filt['mem.util.used'].max()
    mu_min = df_filt['mem.util.used'].min()
    if mu_max != mu_min:
        return (mu_max - mu_min) / 1000
    else:
        return (df_mem[start_index:end_index+2]['mem.util.used'].max() - df_mem[start_index:end_index+2]['mem.util.used'].min()) / 1000


def update_peak_memory_usages(csv, test='optimizations'):
    df = pd.read_csv(csv)

    if (test == 'optimizations'):
        crop_pmus, compute_pmus, mosaic_pmus = [], [], []     
        for index, row in df.iterrows():
            # Get total runtimes of each component from row
            crop_time = int(row['crop_time'])
            compute_time = int(row['compute_time'])
            mosaic_time = int(row['mosaic_time'])

            mem_log_file = f"mem_logs/{row['method']}_{row['library']}_{row['parameter']}_{ast.literal_eval(row['tile_size'])[0]}_{row['run']}.csv"
            crop_pmus.append(get_peak_memory_usage(mem_log_file, 0, crop_time))
            compute_pmus.append(get_peak_memory_usage(mem_log_file, crop_time, compute_time-crop_time))
            mosaic_pmus.append(get_peak_memory_usage(mem_log_file, compute_time, mosaic_time-compute_time-crop_time))

        # Add new columns
        df['crop_pmu'] = crop_pmus
        df['compute_pmu'] = compute_pmus
        df['mosaic_pmu'] = mosaic_pmus

        # Update csv
        df.to_csv(csv, index=False)
    elif (test == 'tile_sizes'):
        pmus = []
        for index, row in df.iterrows():
            # Calculate peak memory usage for each row
            mem_log_file = f"mem_logs/{row['method']}_{row['parameter']}_{ast.literal_eval(row['tile_size'])[0]}_{row['run']}.csv"
            df_mem = pd.read_csv(mem_log_file)
            pmu = (df_mem['mem.util.used'].max() - df_mem['mem.util.used'].min()) / 1000
            pmus.append(pmu)

        # Add new column
        df['peak_memory_usage'] = pmus

        # Update csv
        df.to_csv(csv, index=False)
    elif (test == 'region_changes'):
        pmus = []
        for index, row in df.iterrows():
            # Calculate peak memory usage for each row
            mem_log_file = f"mem_logs/{row['region']}_{row['method']}_{row['parameter']}_{ast.literal_eval(row['tile_size'])[0]}_{row['run']}.csv"
            df_mem = pd.read_csv(mem_log_file)
            pmu = (df_mem['mem.util.used'].max() - df_mem['mem.util.used'].min()) / 1000
            pmus.append(pmu)

        # Add new column
        df['peak_memory_usage'] = pmus

        # Update csv
        df.to_csv(csv, index=False)


def average_together_results(csv, test='optimizations'):
    df = pd.read_csv(csv)

    # Average together with appropiate groupings and save to a new CSV
    if (test == 'optimizations'):
        df_averaged = df.groupby(['method','tile_size','library','parameter'], as_index=False)[['crop_time','compute_time','mosaic_time','crop_pmu','compute_pmu','mosaic_pmu']].mean()
        df_averaged.to_csv(f"averaged_{csv}", index=False)
    elif (test == 'tile_sizes'):
        df_averaged = df.groupby(['tile_size','method','parameter'], as_index=False)[['compute_time','peak_memory_usage']].mean()
        df_averaged.to_csv(f"averaged_{csv}", index=False)
    elif (test == 'process_counts'):
        df_averaged = df.groupby(['method','tile_size','parameter','processes'], as_index=False)[['compute_time']].mean()
        df_averaged.to_csv(f"averaged_{csv}", index=False)
    elif (test == 'region_changes'):
        df_averaged = df.groupby(['region','method','tile_size','parameter'], as_index=False)[['compute_time','peak_memory_usage']].mean()
        df_averaged.to_csv(f"averaged_{csv}", index=False)

###############################
### VISUALIZATION FUNCTIONS ###
###############################

def print_optimization_results(csv):
    df = pd.read_csv(csv)    
    print(tabulate(df.values.tolist(), headers=list(df.columns.values)))


def plot_tile_size_results(csv, terrain_parameter, ylims1, ylims2, zoom_ylims, use_legend=False):
    df = pd.read_csv(csv)
    df_p = df[(df['parameter'] == terrain_parameter)]
    
    tile_sizes = df_p['tile_size'].unique()
    methods = df_p['method'].unique()
    
    execution_times = {}
    pmem_usages = {}
    
    for method in methods:
        execution_times.update({method: df_p[(df_p['method'] == method)]['compute_time'].tolist()})
        pmem_usages.update({method: [val / 1000 for val in df_p[(df_p['method'] == method)]['peak_memory_usage'].tolist()]})
    
    # Configure some plot features
    x = np.arange(len(tile_sizes))  
    width = 0.2
    
    fig, ax = plt.subplots(layout='constrained', figsize=(10, 4))
    ax2 = ax.twinx()
    
    # Plot compute times
    bar_index = 0
    multiplier = 0
    colors = ["sandybrown", "lightblue", "lightgreen", "khaki"]
    hatches = [None,'/','\\','|']
    for attribute, measurement in execution_times.items():
        offset = width * multiplier - 0.1
        rects = ax.bar(x + offset, measurement, width, color=colors[bar_index], hatch=hatches[bar_index], label=attribute)
        multiplier += 1
        bar_index += 1
    
    # Plot peak memory usages
    line_index = 0
    colors = ["red", "blue", "green", "purple"]
    markers = ['*','v','o','^']
    for attribute, measurement in pmem_usages.items():
        ax2.plot(x + width, measurement, linewidth=3, color=colors[line_index], marker=markers[line_index], markersize=10, label=attribute)
        line_index += 1
        
    # Create the zoomed-in axes
    axins = inset_axes(ax, width="30%", height="30%", bbox_to_anchor=(0.04, 0, 1, 0.98), bbox_transform=ax.transAxes, loc='upper left') 
    
    bar_index = 0
    multiplier = 0
    colors = ["sandybrown", "lightblue", "lightgreen", "khaki"]
    for attribute, measurement in execution_times.items():
        offset = width * multiplier - 0.1
        rects = axins.bar(x + offset, measurement, width, color=colors[bar_index], hatch=hatches[bar_index], label=attribute)
        multiplier += 1
        bar_index += 1

    axins.set_xticks(x + width, tile_sizes)
    axins.set_xlim(-0.25, 0.65)
    axins.set_ylim(zoom_ylims[0], zoom_ylims[1])

    # Draw a rectangle to show the zoom area
    mark_inset(ax, axins, loc1=4, loc2=3, fc="none", ec="0.5")
        
    ax.set_xticks(x + width, tile_sizes)
    ax.set_xlabel("Tile Size ((x,y) pixels)", fontweight='bold')
    ax.set_ylabel("Ex. Time (s)", fontweight='bold')
    ax2.set_ylabel("Peak Mem. (GB)", fontweight='bold')
    ax.set_ylim(ylims1[0],ylims1[1])
    ax2.set_ylim(ylims2[0],ylims2[1])
    if use_legend:
        ax.legend(loc='center left', bbox_to_anchor=(0, 1.07), title_fontproperties={'weight':'bold'}, handletextpad=0.4, columnspacing=0.8, ncol=len(methods), frameon=False, title="Execution Time")
        ax2.legend(loc='center left', bbox_to_anchor=(0.515, 1.07), title_fontproperties={'weight':'bold'}, handletextpad=0.4, columnspacing=0.8, ncol=len(methods), frameon=False, title="Peak Memory Usage")

    plt.show()
    
    # Save figure
    plt.savefig(f"imgs/experiment_1_{terrain_parameter}.png")
    plt.clf()
    plt.close()


def plot_process_count_results(csv, terrain_parameter, tile_size, ylims, use_legend=False):
    df = pd.read_csv(csv)
    df_p = df[(df['parameter'] == terrain_parameter) & (df['tile_size'] == f"({tile_size},{tile_size})")]
    
    methods = df_p['method'].unique()
    processes = df_p['processes'].unique()
    
    execution_times = {}
    for method in methods:
        execution_times.update({method: df_p[(df_p['method'] == method)]['compute_time'].tolist()})

    # Configure some plot features
    x = np.arange(len(processes))  
        
    # Plot data
    fig, ax = plt.subplots(layout='constrained', figsize=(8, 4))
    
    # Plot execution times
    line_index = 0
    colors = ["red", "blue"]
    markers = ['*','o']
    for attribute, measurement in execution_times.items():
        ax.plot(x, measurement, linewidth=3, color=colors[line_index], marker=markers[line_index], markersize=10, label=attribute)
        line_index += 1

    ax.set_xticks(x, processes)
    ax.set_xlabel("Number of Processes", fontweight='bold')
    ax.set_ylabel("Ex. Time (s)", fontweight='bold')
    ax.set_ylim(ylims[0],ylims[1])
    if use_legend:
        ax.legend(loc='center left', bbox_to_anchor=(0.32, 1.04), handletextpad=0.4, columnspacing=0.8, ncol=len(methods), frameon=False)

    plt.show()
    
    # Save figure
    plt.savefig(f"imgs/experiment_2_{terrain_parameter}_{tile_size}.png")
    plt.clf()
    plt.close()


def plot_region_change_results(csv, reg, tile_size, ylims1, ylims2, use_legend=False):
    df = pd.read_csv(csv)
    df_p = df[(df['region'] == reg) & (df['tile_size'] == f"({tile_size},{tile_size})")]
    
    methods = df['method'].unique()
    parameters = df['parameter'].unique()
    
    execution_times = {}
    pmem_usages = {}
    
    for method in methods:
        times, pmems = [], []
        for group in PARAMETER_GROUPINGS:
            for param in PARAMETER_GROUPINGS[group]:
                if param in parameters:
                    times.append(df_p[(df_p['method'] == method) & (df_p['parameter'] == param)]['compute_time'].iloc[0])
                    pmems.append(df_p[(df_p['method'] == method) & (df_p['parameter'] == param)]['peak_memory_usage'].iloc[0])
        execution_times.update({method: [float('nan') if val == 0 else val for val in times]})
        pmem_usages.update({method: [float('nan') if val == 0 else val for val in pmems]})
        
    # Configure some plot features
    x = np.arange(len(parameters))  
    width = (1/3)
        
    # Plot data
    fig, ax = plt.subplots(layout='constrained', figsize=(10, 4))
    ax2 = ax.twinx()
    
    # Plot compute times
    bar_index = 0
    multiplier = 0
    colors = ["sandybrown", "lightblue", "lightgreen", "khaki"]
    hatches = [None,'/','\\','|']
    for attribute, measurement in execution_times.items():
        offset = width * multiplier + (1/6)
        rects = ax.bar(x + offset, measurement, width, color=colors[bar_index], hatch=hatches[bar_index], label=attribute)
        multiplier += 1
        bar_index += 1
    
    # Plot peak memory usages
    line_index = 0
    colors = ["red", "blue", "green", "purple"]
    markers = ['*','v','o','^']
    for attribute, measurement in pmem_usages.items():
        ax2.plot(x + width, measurement, linewidth=3, color=colors[line_index], marker=markers[line_index], markersize=10, label=attribute)
        line_index += 1

    # Get abbreviations of terrain parameters for x tick labels
    abbreviated_params, dash_indices = [], []
    dash_index = -1
    for group in PARAMETER_GROUPINGS:
        for param in PARAMETER_GROUPINGS[group]:
            if param in parameters:
                abbreviated_params.append(PARAMETER_CODES[param])
                dash_index += 1
        if (dash_index not in dash_indices) and (dash_index != -1) and ((dash_index != len(parameters)-1)):
            dash_indices.append(dash_index)

    # Add dashed lines between different groupings
    for index in dash_indices:
        x_position = index + 0.8375  # Position between bars
        plt.axvline(x=x_position, color='black', linestyle='--')
    
    ax.set_xticks(x + width, abbreviated_params)
    ax.set_xlabel("Terrain Parameter", fontweight='bold')
    ax.set_ylabel("Ex. Time (s)", fontweight='bold')
    ax2.set_ylabel("Peak Mem. (GB)", fontweight='bold')
    ax.set_ylim(ylims1[0],ylims1[1])
    ax2.set_ylim(ylims2[0],ylims2[1])
    if use_legend:
        ax.legend(loc='center left', bbox_to_anchor=(0.2, 1.07), title_fontproperties={'weight':'bold'}, handletextpad=0.4, columnspacing=0.8, ncol=len(methods), frameon=False, title="Execution Time")
        ax2.legend(loc='center left', bbox_to_anchor=(0.5, 1.07), title_fontproperties={'weight':'bold'}, handletextpad=0.4, columnspacing=0.8, ncol=len(methods), frameon=False, title="Peak Memory Usage")

    plt.show()
    
    # Save figure
    plt.savefig(f"imgs/experiment_3_{reg}_{tile_size}.png")
    plt.clf()
    plt.close()


def print_memory_over_time_results(csv, xlims, ylims):
    df = pd.read_csv(csv)
    
    # Get memory values
    x = np.arange(len(df['mem.util.used']))
    y = [(val/10000) for val in df['mem.util.used'].tolist()]
    
    # Smooth curve
    f = interp1d(x, y, kind='cubic')
    x_smooth = np.linspace(min(x), max(x), 200)
    y_smooth = f(x_smooth)
    
    fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))

    ax.plot(x_smooth,y_smooth)
    ax.set_xlabel("Time (s)", fontweight='bold')
    ax.set_ylabel("Memory Usage (GB)", fontweight='bold')
    ax.set_ylim(ylims[0],ylims[1])
    ax.set_xlim(xlims[0],xlims[1])

    plt.show()
    
    # Save figure
    plt.savefig(f"{csv.replace('mem_logs/','imgs/').replace('.csv','')}_mem_use.png")
    plt.clf()
    plt.close()