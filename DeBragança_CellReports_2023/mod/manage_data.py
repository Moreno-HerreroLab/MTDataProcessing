import numpy as np
import pandas as pd
import mod.manage_pickles as pic


def align_and_classify_data(filepath, del_beads, align_to_step=1, cycles_start_on_step=2):
    """ """
    data, beads_labels = import_data(filepath)
    total = len(beads_labels)
    data, beads_labels = delete_beads(data, beads_labels, del_beads)
    n = len(beads_labels)
    data_aligned, values_align_to = align_data(data, beads_labels, align_to_step)
    # classify the data with respect to the applied F
    data_cycles, data_low_force, data_high_force = classify_data(data_aligned, cycles_start_on_step)

    data_processed = {'aligned': data_aligned,
                      'raw': data_cycles,
                      'high_force': data_high_force,
                      #'low_force': data_low_force,
                      'beads_labels': beads_labels,
                      'z_aligned_to': values_align_to,
                      'total_sample': total,
                      'valid_sample, N': n}

    return data_processed


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
    """
    classifies the data with respect to the experimental part (CYCLES) and the force ('high_force' or 'low_force')
    """
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


def manage_data_from_pickles(folder_path):
    """import and regroup data from pickle"""
    pickledata = pic.import_data_from_pickle(folder_path)
    dataset = filter_according_to_force(pickledata)
    return dataset


def filter_according_to_force(pickledata):
    print('Imported:', end=' ')
    dataset = {}
    conditions = pickledata.keys()
    for cond in conditions:
        print(cond, end=', ')
        filtered, grouped = {}, {}
        merged_by_areas, merged_tot = {}, []
        rep, i = 0, 0
        areas = pickledata[cond].keys()
        for a in areas:
            rep += 1
            areadata = pickledata[cond][a]
            filtered_area, merged = {}, []
            labels = [col for col in areadata['high_force'].columns if 'Z' in col and 'Bead' in col]
            for label in labels:
                data = areadata['high_force'][label]
                filtered_area[f'Bead#{i}'] = data
                merged_tot.append(data)
                merged.append(data)
                i += 1
            filtered[f'rep#{rep}'] = filtered_area
            grouped.update(filtered_area)
            merged_by_areas[f'merged#{rep}'] = pd.concat(merged)    
        filtered['Step'] = areadata['high_force']['Step']
        filtered['Total#Molecules'] = i
        filtered['Time (sg)'] = areadata['high_force']['Time (sg)']
        filtered['Grouped'] = grouped
        filtered['Merged_by_areas'] = merged_by_areas
        filtered['Merged_total'] = pd.concat(merged_tot) 
        dataset[cond] = filtered
    print('done.')
    return dataset

