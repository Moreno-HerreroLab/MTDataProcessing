import seaborn as sns
import numpy as np
import mod.fitting as fit
import mod.process_data as prc
import pandas as pd


def plot_times(axbig, data_to_plot, values_of_interest):
    merge, means = [], []
    for line in data_to_plot:
        cond, _ = line
        tau = np.array(values_of_interest[cond]['time_dwell_merged'])/120
        means += [tau.mean()]
        merge += [tau]
        
    axbig = sns.violinplot(data=merge, inner=None, cut=0, color='grey', scale='width', scale_hue=True)
    for violin in axbig.collections:
        violin.set_linewidth(1)
        violin.set_alpha(0.7)
        violin.set_edgecolor('grey')
    
    i=0
    for mean in means:
        axbig.scatter(x=i, y=mean, marker='o', s=25, color='black', edgecolor='white')
        axbig.annotate(f'{round(mean,2)} s', xycoords='data', xy=(i, mean+0.5), size=10, ha='center', va='center')
        i+=1
    axbig.set(ylim=(None, None), ylabel='Dwell time, t (s)')


def plot_distribution(ax, line, dataset, bw, sigmas, beta, area=None, window=[False, None]):
    # data
    cond, label = line
    w, window = window[0], window[1]
    if area != None:
        data = dataset[cond]['Merged_by_areas'][f'merged#{area}']
        molecules = len(dataset[cond][f'rep#{area}'])
    elif w:
        data = []
        #print(window)
        for bead in window: data.append(dataset[cond]['Grouped'][bead])
        data = pd.concat(data)
        molecules = len(window)
    else: 
        data = dataset[cond]['Merged_total']
        molecules = len(dataset[cond]['Grouped'])
    # get histogram parameters for data    
    #bins = np.arange(min(data), max(data), bw)
    bins = np.arange(-0.65, 0.15, bw)
    counts, bins_edges = np.histogram(data, bins=bins, density=True)
    bins_center = np.array(prc.get_bin_central_points(bins_edges))
    # plot entire histogram 
    ax.hist(bins_center, bins, weights=counts, fc='tomato', ec='white', alpha=0.5, lw=1)
    # gaussian fitting of the individual expected peaks
    mask_extended = (-0.05 <= bins_center) 
    mask_synaptic = np.logical_and(-0.3 <= bins_center, bins_center <= -0.1)
    common = ax, bw
    ext_percent = fit.fit_peak(bins_center, counts, (mask_extended, 0-beta, 0+beta), common, sigmas) # Extended peak 
    syn_percent = fit.fit_peak(bins_center, counts, (mask_synaptic, -0.2-beta, -0.2+beta), common, sigmas) # synaptic peak
    percentages = ext_percent, syn_percent
    # annotations ---
    bbox_props = dict(boxstyle='round', fc='white', ec='0.1', alpha=0.9)
    ax.text(x=-0.63, y=3.1, s=f'{label}', ha='left', va='center', size=10, bbox=bbox_props)
    if w: ax.annotate(f'{window[0:3]}[...]', xy=(-0.63, 0.8), size=6, ha='left', va='center' )
    ax.annotate(f'N = {molecules}', xy=(-0.63, 0.5), size=10, ha='left', va='center')
    ax.set(ylabel='Density', axisbelow=True)
    ax.grid(color='gray', alpha=0.5, linestyle='--')
    return percentages