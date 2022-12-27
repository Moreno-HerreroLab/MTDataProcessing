import pandas as pd
import numpy as np


def extracting_values_of_interest(dataset, threshold_um, roll_window, factor):
    conditions, i = dataset.keys(), 1
    inputs = threshold_um, roll_window, factor
    extracted_values = {'Inputs':   {'threshold_um': threshold_um, 
                                    'roll_window': roll_window,
                                    'points_for_means': roll_window*factor}}
    for cond in conditions:
        print(f'Working on condition {i}/{len(conditions)} - {cond}...', end='')
        extracted_values[cond] = extract_within_condition(dataset[cond], inputs)
        i+=1
        print(' done.')
    print('Finished!')
    return extracted_values


def extract_within_condition(data, inputs):
    data_all = sort_all_beads(data)
    step_values = np.unique(data['Step'])
    extracted_dict = {}
    extracted_values = create_lists_and_counts()
    for bead in data_all.keys():
        extracted_dict[bead] = {}
        deltas, dwell, num_events = create_lists_and_counts()
        for cycle in range(0, len(step_values)):
            # filter data per cycle
            mask = (data['Step'] == step_values[cycle])
            points = data_all[bead][mask]
            # get values
            down = downsample_data(points, inputs)
            extracted_dict[bead][f'Cycle#{cycle}'] = {'data': points, 'downsampled': down}
            dwell, num_events = extract_rupture_values(points, inputs, dwell, num_events)
            deltas = extract_delta_ini_and_fin(points, inputs, deltas)
        store_extracted_values_for_bead(extracted_dict[bead], deltas, dwell, num_events)
        extracted_values = merge_extracted_values_within_condition(extracted_values, deltas, dwell, num_events)
    store_merged_extracted_values(extracted_dict, extracted_values)
    return extracted_dict


def sort_all_beads(data):
    repetitions = [i for i in data.keys() if 'rep' in i]
    data_all = {}
    for rep in repetitions:
        for bead in data[rep].keys():
            data_all[bead] = data[rep][bead]
    return data_all


def create_lists_and_counts():
    # deltas
    deltas = [[], []]
    # rupture delta and time
    dwell = [[], []]
    # events
    num_events = 0
    return deltas, dwell, num_events


def extract_rupture_values(points, inputs, dwell, num_events):
    threshold_um, roll_window, factor = inputs
    points_for_mean = factor*roll_window
    delta_dwell, time_dwell = dwell
    # --
    ini_frame = points.index[0]
    downsampled = downsample_data(points, inputs)
    # get the time of rupture
    mask = downsampled >= -threshold_um
    if len(downsampled[mask]) == 0:
        fin_frame = points.index[-1]
    else:
        fin_frame = downsampled[mask].index[0]
    tau_frame = fin_frame - ini_frame
    time_dwell += [tau_frame]
    # get the delta at the rupture time
    if tau_frame > roll_window:
        num_events += 1
        if tau_frame > points_for_mean:
            ini_rupture = tau_frame - points_for_mean
        else: ini_rupture = 0
        delta = points[ini_rupture:tau_frame].mean()
    else:
        delta = points[:points_for_mean].mean()
    delta_dwell += [delta]
    return dwell, num_events


def downsample_data(points, inputs):
    _, roll_window, _ = inputs
    df = pd.DataFrame(points)
    downsampled = df.rolling(window=roll_window, min_periods=1, center=True).mean()
    downsampled = downsampled[downsampled.columns[0]]
    return downsampled


def extract_delta_ini_and_fin(points, inputs, delta_values):
    """get the delta at the beginning of the cycle"""
    delta_ini, delta_fin = delta_values
    _, roll_window, factor = inputs
    points_for_mean = factor*roll_window
    # --
    ini = points[:points_for_mean].mean()
    fin = points[-points_for_mean:].mean()
    delta_ini += [ini]
    delta_fin += [fin]
    return delta_values


def store_extracted_values_for_bead(d, deltas, dwell, num_events):
    """ save the extracted values to the dictionary """
    delta_ini, delta_fin = deltas
    delta_dwell, time_dwell = dwell
    # --
    d['delta_ini'] = delta_ini
    d['delta_fin'] = delta_fin
    d['delta_dwell'] = delta_dwell
    d['time_dwell'] = time_dwell
    d['events'] = [num_events]


def merge_extracted_values_within_condition(extracted_values, deltas, dwell, num_events):
    """ save the extracted values to the merged datasets"""
    deltas_merged, rupture_merged, molecules_with_events = extracted_values
    # ---
    delta_ini, delta_fin = deltas
    delta_dwell, time_dwell = dwell
    # -- 
    deltas_merged[0] += delta_ini  # delta ini list
    deltas_merged[1] += delta_fin  # delta fin list
    # --
    rupture_merged[0] += delta_dwell   # dwell delta list
    rupture_merged[1] += time_dwell    # dwell time list
    # --
    if num_events > 0: molecules_with_events += 1
    # --
    extracted_values = deltas_merged, rupture_merged, molecules_with_events
    return extracted_values


def store_merged_extracted_values(d, extracted_values):
    """ save the extracted values to the dictionary """
    deltas_merged, rupture_merged, molecules_with_events = extracted_values
    delta_ini_merged, delta_fin_merged = deltas_merged
    delta_rupture_merged, time_rupture_merged = rupture_merged
    # --
    d['delta_ini_merged'] = delta_ini_merged
    d['delta_fin_merged'] = delta_fin_merged
    d['delta_dwell_merged'] = delta_rupture_merged
    d['time_dwell_merged'] = time_rupture_merged
    d['molecules_with_events'] = molecules_with_events
