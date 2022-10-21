import matplotlib.pyplot as plt
import extract_values
import pandas as pd
import numpy as np
import pickle
import os


def def_fig_parameters():
    plt.rc('axes', titlesize=14, labelsize=12)  # font size of the axes title and of the x and y labels
    plt.rc('xtick', labelsize=10, direction='in', top=True)  # font size of the tick labels
    plt.rc('ytick', labelsize=10, direction='in', right=True)  # font size of the tick labels
    plt.rc('legend', fontsize=10, markerscale=10, loc='lower left')  # legend font size
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    
def align_classify_savepickle(folder, condition, results, save=False):
    exp = {}
    area = 0
    for line in results:
        R, del_beads = line
        area += 1
        path = os.path.join(folder, 'data_corrected', f'results_{R}_PAR_corr.dat')
        data, beads_labels = import_data(path)
        data, beads_labels = delete_beads(data, beads_labels, del_beads)  # delete problematic beads (lost while tracking...)
        data, align_to = align_data(data, beads_labels, align_to_step=1) # sustract the max extension to each point of each bead 
        data_cycles, data_low_force, data_high_force = classify_data(data, cycles_start_on_step=2) # with respect to cycles and applied F
        data_aligned_classified = {'aligned': data,  # all data -> incubation, baseline and cycles
                                   'cycles': data_cycles,  # all data from the force cycles part of the experiment
                                   'high_force': data_high_force,  # data from the higher force steps
                                   'low_force': data_low_force,  # data for the lower force steps
                                   'beads_labels': beads_labels, 
                                   'z_aligned_to': align_to  # maximum extensions for each bead
                                  }
        
        exp[f'area{area}'] = data_aligned_classified
        if save: save_data_to_pickle(condition, folder, area, data_aligned_classified)
    print('Finished!')
    return data_aligned_classified


def import_data(path):
    data = pd.read_csv(path, delimiter='\t', header=0)
    beads_labels = [col for col in data.columns if 'Bead' in col and 'Z' in col]
    return data, beads_labels


def align_data(data, beads_labels, align_to_step):
    data_aligned = data.copy(deep=True)
    aligned_with_respect_to = {}
    for label in beads_labels:
        ext = np.mean(data[label][data['Step'] == align_to_step])
        aligned_with_respect_to[label] = ext
        diff = 0 - ext
        data_aligned[label] = data_aligned[label] + diff
    return data_aligned, aligned_with_respect_to


def classify_data(data, cycles_start_on_step):
    mask = data['Step'] >= cycles_start_on_step
    data_cycles = data[mask].copy(deep=True)

    pos1 = data['Shift pos (mm)'][data['Step'] == cycles_start_on_step]
    pos2 = data['Shift pos (mm)'][data['Step'] == cycles_start_on_step + 1]
    position_low_force = round(pos1.values[0], 1) if pos1.values[0] > pos2.values[0] else round(pos2.values[0], 1)

    data_low_force = data_cycles[round(data_cycles['Shift pos (mm)'], 1) == position_low_force]
    data_high_force = data_cycles[round(data_cycles['Shift pos (mm)'], 1) != position_low_force]
    return data_cycles, data_low_force, data_high_force


def delete_beads(data, beads_labels, del_beads):
    for number in del_beads:
        columns = [f'Bead{number}-X(um)', f'Bead{number}-Y(um)', f'Bead{number}-Z(um)']
        data.drop(columns, axis=1, inplace=True)
        beads_labels.remove(columns[-1])
    return data, beads_labels


def regroup_data_within_condition(dataset):
    """regrouping the data for a particular condition from different dates into a single dataset"""
    dataset_regrouped = {}
    for condition in list(dataset.keys()):
        print(condition)
        data_regrouped = {}
        i = 0
        dates = list(dataset[condition].keys())
        for date in dates:
            data_unmerged = dataset[condition][date]['HF']
            columns = list(data_unmerged.keys())
            labels = [col for col in columns if 'Z' in col and 'Bead' in col]
            print(len(labels), 'beads', date)
            for label in labels:
                data_regrouped[f'Bead#{i}'] = dataset[condition][date]['HF'][label]
                i += 1
        data_regrouped['Step'] = dataset[condition][dates[0]]['HF']['Step']
        data_regrouped['Time (sg)'] = dataset[condition][dates[0]]['HF']['Time (sg)']
        print('total:', len(data_regrouped.keys()) - 2, 'beads')
        dataset_regrouped[condition] = data_regrouped
    return dataset_regrouped



def save_data_to_pickle(condition, date, area, data_aligned):
    if not os.path.isdir('Pickles'): os.mkdir('Pickles')
    path = os.path.join('Pickles', condition)
    data_to_file = {f'{date}_area{area}': {'cycles': data_aligned['cycles'],
                                           'HF': data_aligned['high_force'],
                                           'LF': data_aligned['low_force']}}
    if os.path.exists(path + '.pickle'):
        # Load data in dictionary
        with open(path + '.pickle', 'rb') as file_handler:
            saved_data = pickle.load(file_handler)
        print(f'Loading saved dictionary.')
        print(f'Content: {list(saved_data.keys())}')
        # save new data to the dictionary
        saved_data[f'{date}_area{area}'] = data_to_file[f'{date}_area{area}']
        data_to_file = saved_data
    else:
        print('Creating new dictionary.')
    with open(path + '.pickle', 'wb') as file_handler:
        pickle.dump(data_to_file, file_handler)
    file_handler.close()
    print(list(data_to_file.keys()))

 
def process_and_merge(dataset, conditions, threshold_um=0.05, roll_window=10, factor=3):
    dataset_processed, dataset_merged, i = {}, {}, 0
    for condition in conditions:
        i+=1
        print(f'Condition {i}/{len(conditions)} - {condition}...', end='')
        # data processing
        data = dataset[condition]
        data_processed = extract_values_of_interest(data, threshold_um, roll_window, factor)
        dataset_processed[condition] = data_processed
        # data merging
        num_of_cycles = dataset_processed[condition]['Info']['#Cycles']
        beads_labels = list(dataset_processed[condition].keys())
        beads_labels = [label for label in beads_labels if 'Bead' in label]
        data = []
        for label in beads_labels:
            for cycle in range(0, num_of_cycles):
                bead_data = dataset_processed[condition][label][f'Cycle#{cycle}']['data']
                data.append(bead_data)
        data = pd.concat(data)
        dataset_merged[condition] = data
        # finish condition
        print(' Done.')
    print('Finished!')    
    return dataset_processed, dataset_merged
    

def manage_data_from_pickles(folder_path):
    """import and regroup data from pickle"""
    dataset = import_data_from_pickle(folder_path)
    dataset_regrouped = regroup_data_within_condition(dataset)
    return dataset_regrouped

    
def import_data_from_pickle(folder_path):
    dataset = {}
    for filename in os.listdir(folder_path):
        if filename.lower().endswith('.pickle'):
            condition = filename.split('.')[0]
            filepath = os.path.join(folder_path, filename)
            with open(filepath, 'rb') as file_handler:
                file_data = pickle.load(file_handler)
            dataset[condition] = file_data
        else:
            pass
    return dataset

def regroup_data_within_condition(dataset):
    """regrouping the data for a particular condition from different dates into a single dataset"""
    dataset_regrouped = {}
    for condition in list(dataset.keys()):
        print(condition)
        data_regrouped = {}
        i = 0
        dates = list(dataset[condition].keys())
        for date in dates:
            data_unmerged = dataset[condition][date]['HF']
            columns = list(data_unmerged.keys())
            labels = [col for col in columns if 'Z' in col and 'Bead' in col]
            print(len(labels), 'beads', date)
            for label in labels:
                data_regrouped[f'Bead#{i}'] = dataset[condition][date]['HF'][label]
                i += 1
        data_regrouped['Step'] = dataset[condition][dates[0]]['HF']['Step']
        data_regrouped['Time (sg)'] = dataset[condition][dates[0]]['HF']['Time (sg)']
        print('total:', len(data_regrouped.keys()) - 2, 'beads')
        dataset_regrouped[condition] = data_regrouped
    return dataset_regrouped



def merge_data(data_low_force, data_high_force, beads_labels):
    """ """
    global_low_force, global_high_force = [], []
    for label in beads_labels:
        global_low_force.append(data_low_force[label])
        global_high_force.append(data_high_force[label])

    global_low_force = pd.concat(global_low_force)
    global_high_force = pd.concat(global_high_force)
    return global_low_force, global_high_force


def extract_values_of_interest(data, threshold_um=0.05, roll_window=10, factor=3):
    data_processed = extract_values.extract_times_and_deltas(data, threshold_um, roll_window, factor)
    return data_processed


def get_bin_central_points(bin_edges):
    bin_center = []
    idx_max = len(bin_edges[:-1])
    for i in range(0, idx_max):
        mean = (bin_edges[i + 1] + bin_edges[i]) / 2
        bin_center += [mean]
    return bin_center

def gaussian_fit(x_to_fit, y_to_fit, xrange, expected, x_bins, bw, axs):
    from pylab import *  # for the gaussian curve fitting
    from scipy.optimize import curve_fit  # for the gaussian curve fitting
    def gauss(x, mu, sigma, A):
        return A*exp(-(x-mu)**2/(2*sigma**2))
    params, cov = curve_fit(gauss, x_to_fit, y_to_fit, expected)
    # plot the gaussian function
    x, y = xrange, gauss(xrange,*params)
    axs[0].plot(x, y, 'k-', linewidth=1)
    ymax = y.max()
    xmax = xrange[y==ymax][0]
    axs[0].annotate(f'{round(xmax, 4)} µm', xy=(xmax, 1),
                      xytext=(xmax, 3), size=10, ha='center', va='center', rotation=90)
                      #arrowprops=dict(arrowstyle='->',lw=1))
    print(f'mu: {params[0]}, sigma: {params[1]}, A: {params[2]}')
    # find the area under the curve
    x, y = x_bins, gauss(x_bins,*params)
    #axs[0].scatter(x, y, s=1, color='k')
    percent = int(sum(y*bw)*100)
    print('area_under_curve=', percent)
    axs[0].annotate(f'{percent} %', xy=(params[0], 0.5), xycoords='data', size=10, ha='center', va='center')

    
def find_maxima_hist_and_annotate(counts, bin_center, round_to=4):
    ymax = counts.max()
    xmax = bin_center[counts==ymax][0]
    return xmax, ymax

def plot_distribution(dataset_merged, dataset_processed, line, coords, axs):
    condition, label = line
    data = dataset_merged[condition]
    molecules = dataset_processed[condition]['Info']['#Molecules']

    bw = 0.015 
    bins = np.arange(min(data), max(data), bw)
    counts, bin_edges = np.histogram(data, bins=bins, density=True)
    bin_center = np.array(get_bin_central_points(bin_edges))
    
    mask_extended = -0.05 <= bin_center 
    mask_synaptic = np.logical_and(-0.25 <= bin_center, bin_center <= -0.15)
    mask_bridged = bin_center < -0.05
    counts_extended, bin_center_extended = counts[mask_extended], bin_center[mask_extended]
    counts_bridged, bin_center_bridged = counts[mask_bridged], bin_center[mask_bridged]
    counts_synaptic, bin_center_synaptic = counts[mask_synaptic], bin_center[mask_synaptic]
    
    axs[0].hist(bin_center_extended, bins, weights=counts_extended,
                color='tomato', edgecolor='white', alpha=0.5, lw=1, label='extended DNA')
    axs[0].hist(bin_center_bridged, bins, weights=counts_bridged,
                color='tomato', edgecolor='white', alpha=0.5, lw=1, label='bridged DNA')
    try:
        # Common parameters
        sigma = .05
        xrange = np.linspace(-0.6, 0.05, 5000)

        # extended molecule peak
        xmax_e, ymax_e = find_maxima_hist_and_annotate(counts=counts_extended, bin_center=bin_center_extended)
        gaussian_fit(x_to_fit=bin_center_extended,y_to_fit=counts_extended, xrange=xrange, 
                     expected=(xmax_e, sigma, ymax_e), x_bins=bin_center, bw=bw, axs=axs)
        
        # synapsis Peak
        xmax_s, ymax_s = find_maxima_hist_and_annotate(counts=counts_synaptic, bin_center=bin_center_synaptic)
        gaussian_fit(x_to_fit=bin_center_synaptic, y_to_fit=counts_synaptic, xrange=xrange, 
                     expected=(xmax_s, sigma, ymax_s), x_bins=bin_center, bw=bw, axs=axs)
        
    except: 
        print('unable to fit')
            
    # annotations ---
    bbox_props = dict(boxstyle='round', fc='white', ec='0.1', alpha=0.9)
    axs[0].text(x=-0.63, y=3.2, s=f'{label}', ha='left', va='center', size=10, bbox=bbox_props)
    axs[0].annotate(f'N = {molecules}', xy=(-0.63, 0.5), size=10, ha='left', va='center')
    axs[0].set(ylabel='Density', axisbelow=True)
    axs[0].grid(color='gray', alpha=0.5, linestyle='--')
    
def get_contribution_percentage(all_counts, bw, axs):
    counts, counts_extended, counts_bridged = all_counts
    area=sum(counts*bw)
    distribution_all = int(area*100)
    distribution_extended = int(sum(counts_extended*bw)*100)
    distribution_bridged = int(sum(counts_bridged*bw)*100)
    print(f'Condition: Norm (%) =',distribution_all, ' Extended (%) =', distribution_extended, ' Not Extended (%) =', distribution_bridged)

    axs[0].annotate(f'{distribution_bridged} %', xy=(-0.35, 0.5), xycoords='data', size=10, ha='center', va='center', color='red')
    axs[0].annotate(f'{distribution_extended} %', xy=(0.1, 0.5), xycoords='data', size=10, ha='center', va='center', color='red')