from matplotlib import pyplot as plt
import pandas as pd
import numpy as np 

def visualize_tile_compute(csv_file, variable, plot_title):
    # Read CSV into pandas dataframe
    data = pd.read_csv(csv_file)
    df = pd.DataFrame(data)

    # Variable options are average_compute_time_per_tile and average_compute_mem_usage_per_tile
    tile_count = df['tile_count']

    slope = df[df['parameter'] == 'slope'][variable]
    aspect = df[df['parameter'] == 'aspect'][variable]
    hillshade = df[df['parameter'] == 'hillshade'][variable]
    pfc = df[df['parameter'] == 'plan_curvature'][variable]
    plc = df[df['parameter'] == 'profile_curvature'][variable]
    ci = df[df['parameter'] == 'convergence_index'][variable]

    # Convert memory usage to KB
    if variable == 'average_compute_mem_usage_per_tile':
        slope = slope / 1024
        aspect = aspect / 1024
        hillshade = hillshade / 1024
        pfc = pfc / 1024
        plc = plc / 1024
        ci = ci / 1024
    
    # set width of bar 
    barWidth = 0.1
    fig = plt.subplots(figsize=(12, 8)) 
    
    # Set position of bar on X axis 
    br1 = np.arange(len(slope)) 
    br2 = [x + barWidth for x in br1] 
    br3 = [x + barWidth for x in br2] 
    br4 = [x + barWidth for x in br3] 
    br5 = [x + barWidth for x in br4] 
    br6 = [x + barWidth for x in br5] 
    
    # Make the plot
    plt.bar(br1, slope, color ='r', width = barWidth, edgecolor ='grey', label ='Slope') 
    plt.bar(br2, aspect, color ='g', width = barWidth, edgecolor ='grey', label ='Aspect') 
    plt.bar(br3, hillshade, color ='b', width = barWidth, edgecolor ='grey', label ='Hillshade') 
    plt.bar(br4, pfc, color ='y', width = barWidth, edgecolor ='grey', label ='Profile Curv') 
    plt.bar(br5, plc, color ='m', width = barWidth, edgecolor ='grey', label ='Plan Curv') 
    plt.bar(br6, ci, color ='c', width = barWidth, edgecolor ='grey', label ='Conv Index') 
    
    # Adding labels
    ylabel = 'Execution Time (s)' if variable == 'average_compute_time_per_tile' else 'Peak Memory Usage (KB)'
    plt.xlabel('Tile Count', fontweight ='bold', fontsize = 15)
    plt.xticks([r + barWidth for r in range(len(slope))], tile_count.unique())
    plt.ylabel(ylabel, fontweight ='bold', fontsize = 15)
    plt.title(plot_title, fontweight ='bold', fontsize = 15)
    
    plt.legend()
    plt.show() 