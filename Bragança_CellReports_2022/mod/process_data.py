import os
import mod.extract_values as extract
import mod.details as det
import mod.manage_data as mani
from mod.figparams import *
import numpy as np


def process(info, fov, align_to_step, cycles_start_on_step):
    folder, condition, date, flowcell  = info
    details = det.import_experimental_details(folder)
    data = {}
    for line in fov:
        area, R, del_beads = line
        print(area, end=' ')
        # ---
        path = os.path.join(folder, date,'data_corrected', f'results_{R}_PAR_corr.dat')
        data[area] = mani.align_and_classify_data(path, del_beads, align_to_step, cycles_start_on_step)
        # ---
        details = det.save_experimental_details(info, data, line, details)
        # ---
        print('done.')
    print('Finished processing.')
    return data, details


def get_bin_central_points(bins_edges):
    bins_center = []
    idx_max = len(bins_edges[:-1])
    for i in range(0, idx_max):
        mean = (bins_edges[i + 1] + bins_edges[i]) / 2
        bins_center += [mean]
    return bins_center


def find_coord_hist_maxima(bins_center, counts):
    y_peak = counts.max()
    x_peak = bins_center[counts==y_peak][0]
    return x_peak, y_peak


"""

def get_contribution_percentage(all_counts, bw): # NOT BEING USED 
    counts, counts_extended, counts_bridged = all_counts
    area=sum(counts*bw)
    distribution_all = int(area*100)
    distribution_extended = int(sum(counts_extended*bw)*100)
    distribution_bridged = int(sum(counts_bridged*bw)*100)
"""