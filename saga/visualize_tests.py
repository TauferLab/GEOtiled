from matplotlib import pyplot as plt
import pandas as pd
import numpy as np 

def visualize_tile_compute(csv_file, variable, plot_title):
    # Read CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    pre_df = pd.DataFrame(data)
    df = pre_df[pre_df['tile_count'] < 101]

    # Variable options are average_compute_time_per_tile and average_compute_mem_usage_per_tile
    tile_count = df['tile_count']

    aspect = df[df['parameter'] == 'aspect'][variable]
    hillshade = df[df['parameter'] == 'hillshade'][variable]
    pfc = df[df['parameter'] == 'profile_curvature'][variable]
    ci = df[df['parameter'] == 'convergence_index'][variable]

    aspect_std = df[df['parameter'] == 'aspect']['avg_comp_time_std']
    hillshade_std = df[df['parameter'] == 'hillshade']['avg_comp_time_std']
    pfc_std = df[df['parameter'] == 'profile_curvature']['avg_comp_time_std']
    ci_std = df[df['parameter'] == 'convergence_index']['avg_comp_time_std']

    # Convert memory usage to KB
    if variable == 'avg_mem_usage_per_tile':
        mb_conversion = (1024*1024)
        aspect = aspect / mb_conversion
        hillshade = hillshade / mb_conversion
        pfc = pfc / mb_conversion
        ci = ci / mb_conversion

        aspect_std = df[df['parameter'] == 'aspect']['avg_mem_usage_std'] / mb_conversion
        hillshade_std = df[df['parameter'] == 'hillshade']['avg_mem_usage_std'] / mb_conversion
        pfc_std = df[df['parameter'] == 'profile_curvature']['avg_mem_usage_std'] / mb_conversion
        ci_std = df[df['parameter'] == 'convergence_index']['avg_mem_usage_std'] / mb_conversion
    
    # set width of bar 
    barWidth = 0.1
    fig = plt.subplots(figsize=(12, 8)) 
    
    # Set position of bar on X axis 
    br1 = np.arange(len(aspect)) 
    br2 = [x + barWidth for x in br1] 
    br3 = [x + barWidth for x in br2] 
    br4 = [x + barWidth for x in br3] 
    
    # Make the plot 
    plt.bar(br1, aspect, color ='r', width = barWidth, edgecolor ='grey', label ='Aspect') 
    plt.errorbar(br1, aspect, yerr=aspect_std, fmt='k', linestyle='')
    plt.bar(br2, hillshade, color ='g', width = barWidth, edgecolor ='grey', label ='Hillshade')
    plt.errorbar(br2, hillshade, yerr=hillshade_std, fmt='k', linestyle='')
    plt.bar(br3, pfc, color ='b', width = barWidth, edgecolor ='grey', label ='Profile Curv') 
    plt.errorbar(br3, pfc, yerr=pfc_std, fmt='k', linestyle='')
    plt.bar(br4, ci, color ='y', width = barWidth, edgecolor ='grey', label ='Conv Index') 
    plt.errorbar(br4, ci, yerr=ci_std, fmt='k', linestyle='')
    
    # Adding labels
    ylabel = 'Execution Time (s)' if variable == 'avg_comp_time_per_tile' else 'Peak Memory Usage (MB)'
    plt.xlabel('Tile Count', fontweight ='bold', fontsize = 15)
    plt.xticks([r + barWidth for r in range(len(aspect))], df['tile_count'].unique())
    plt.ylabel(ylabel, fontweight ='bold', fontsize = 15)
    plt.title(plot_title, fontweight ='bold', fontsize = 15)
    
    plt.legend()
    plt.show() 

def visualize_crop_test(csv_file, plot_title):
    # Read CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    tc = df['tile_count']
    et = df['execution_time']

    # Figure Size
    fig = plt.figure(figsize =(10, 7))

    # Create bar plot
    barWidth = 0.75
    bar = np.arange(len(tc)) 
    plt.bar(bar, et, color ='c', width = barWidth, edgecolor ='grey', label ='Execution Time') 

    # Labels
    plt.xticks([r for r in range(len(tc))], tc)
    plt.xlabel('Tile Count', fontweight ='bold', fontsize = 15)
    plt.ylabel('Execution Time (s)', fontweight ='bold', fontsize = 15)
    plt.title(plot_title, fontweight ='bold', fontsize = 15)
    
    # Show Plot
    plt.show()

def visualize_timing_test(csv_file, plot_title):
    # Read CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    tc = df['tile_count']
    et = df['execution_time']

    # Figure Size
    fig = plt.figure(figsize =(10, 7))

    # Create bar plot
    barWidth = 0.75
    bar = np.arange(len(tc)) 
    plt.bar(bar, et, color ='c', width = barWidth, edgecolor ='grey', label ='Execution Time (s)') 

    # Labels
    plt.xticks([r for r in range(len(tc))], tc)
    plt.xlabel('Tile Count', fontweight ='bold', fontsize = 15)
    plt.ylabel('Execution Time (s)', fontweight ='bold', fontsize = 15)
    plt.title(plot_title, fontweight ='bold', fontsize = 15)
    
    # Show Plot
    plt.show()

def visualize_memory_test(csv_file, plot_title):
    # Read CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    tc = df['tile_count']
    pmu = df['peak_mem_usage'] / (1024*1024)

    # Figure Size
    fig = plt.figure(figsize =(10, 7))

    # Create bar plot
    barWidth = 0.75
    bar = np.arange(len(tc)) 
    plt.bar(bar, pmu, color ='m', width = barWidth, edgecolor ='grey', label ='Peak Memory Usage (MB)') 

    # Labels
    plt.xticks([r for r in range(len(tc))], tc)
    plt.xlabel('Tile Count', fontweight ='bold', fontsize = 15)
    plt.ylabel('Peak Memory Usage (MB)', fontweight ='bold', fontsize = 15)
    plt.title(plot_title, fontweight ='bold', fontsize = 15)
    
    # Show Plot
    plt.show()