import pandas as pd


def extract_times_and_deltas(data, threshold_um, roll_window, factor):
    """extract the times and deltas of interest"""
    inputs = [threshold_um, roll_window, factor]
    data_processed, steps, labels = store_info(data)
    extracted_merged = create_merged_lists_and_counts()
    for bead in labels:
        data_processed[bead] = {}
        delta_values, rupture_values, events_values = refresh_lists_and_counts()
        for cycle in range(0, len(steps)):
            # filter data per cycle
            step_num = steps[cycle]  # steps is 2*num of cycle
            mask = data['Step'] == step_num
            points = data[bead][mask]
            # store points
            downsampled = downsample_data(points, roll_window)
            data_processed[bead][f'Cycle#{cycle}'] = {'data': points,
                                                      'downsampled': downsampled}
            # get time_rupture and delta_rupture
            rupture_values, events_values = extract_rupture_values(points, inputs, rupture_values, events_values)
            # get deltas
            delta_values = extract_delta_ini_and_fin(points, inputs, delta_values)
        store_bead_extracted_values(data_processed[bead], delta_values, rupture_values, events_values)
        extracted_merged = merge_extracted_values(extracted_merged, delta_values, rupture_values, events_values)
    store_merged_extracted_values(data_processed, extracted_merged)
    return data_processed


def store_info(d):
    labels = get_beads_labels(d)
    indexes, cycles, steps = get_steps_and_indexes(d)
    info = {'Labels': labels,
            '#Molecules': len(labels),
            '#Cycles': cycles,
            'Steps': steps,
            'Indexes': indexes}
    data_processed = {'Info': info}
    return data_processed, steps, labels


def get_beads_labels(data):
    labels = [col for col in data if 'Bead' in col]
    return labels


def get_steps_and_indexes(data):
    """ get the indexes that corresponds to each new step """
    steps = []
    indexes = []
    for element in data['Step']:
        if element not in steps:
            steps = steps + [element]
            mask = data['Step'] == element
            indexes = indexes + [data['Step'][mask].index[0]]
    cycles = len(steps)
    return indexes, cycles, steps


def create_merged_lists_and_counts():
    # create lists for delta values
    delta_ini_merged, delta_fin_merged = [], []
    deltas_merged = [delta_ini_merged, delta_fin_merged]
    # create lists for rupture values
    delta_rupture_merged, time_rupture_merged = [], []
    rupture_merged = [delta_rupture_merged, time_rupture_merged]
    # create variables for events values
    molecules_with_events, unbroken_events = [0], [0]
    events_merged = [molecules_with_events, unbroken_events]
    # joined list
    extracted_merged = [deltas_merged, rupture_merged, events_merged]
    return extracted_merged


def refresh_lists_and_counts():
    # create lists for deltas
    delta_ini, delta_fin = [], []
    delta_values = [delta_ini, delta_fin]
    # create lists for rupture values
    time_rupture, delta_rupture = [], []
    rupture_values = [time_rupture, delta_rupture]
    # lists for events
    events, unbroken = [0], [0]
    events_values = [events, unbroken]
    return delta_values, rupture_values, events_values


def extract_rupture_values(points, inputs, rupture_values, events_values):
    threshold_um, roll_window, factor = inputs
    points_for_mean = factor*roll_window
    time_rupture, delta_rupture = rupture_values
    events, unbroken = events_values
    # --
    ini_frame = points.index[0]
    downsampled = downsample_data(points, roll_window)
    # get the time of rupture
    mask = downsampled >= -threshold_um
    if len(downsampled[mask]) == 0:
        # when it lasts the 10 s, unbroken event
        fin_frame = points.index[-1]
        unbroken[0] += 1
    else:
        fin_frame = downsampled[mask].index[0]
    tau_frame = fin_frame - ini_frame
    time_rupture += [tau_frame]
    # get the delta rupture
    if tau_frame > roll_window:
        # if true, then it is an event
        events[0] += 1
        if tau_frame > points_for_mean:
            ini = tau_frame - points_for_mean
            delta = points[ini:tau_frame].mean()
        else:
            delta = points[:tau_frame].mean()
    else:
        delta = points[:points_for_mean].mean()
    delta_rupture += [delta]
    return rupture_values, events_values


def downsample_data(points, roll_window):
    df = pd.DataFrame(points)
    downsampled = df.rolling(window=roll_window, min_periods=1, center=True).mean()
    downsampled = downsampled[downsampled.columns[0]]
    return downsampled


def extract_delta_ini_and_fin(points, inputs, delta_values):
    """get the delta at the beginning of the cycle"""
    delta_ini, delta_fin = delta_values
    threshold_um, roll_window, factor = inputs
    points_for_mean = factor*roll_window
    # --
    ini = points[:points_for_mean].mean()
    fin = points[-points_for_mean:].mean()
    delta_ini += [ini]
    delta_fin += [fin]
    return delta_values


def store_bead_extracted_values(d, delta_values, rupture_values, events_values):
    """ save the extracted values to the dictionary """
    delta_ini, delta_fin = delta_values
    time_rupture, delta_rupture = rupture_values
    events, unbroken = events_values
    # --
    d['delta_ini'] = delta_ini
    d['delta_fin'] = delta_fin
    d['delta_rupture'] = delta_rupture
    d['time_rupture'] = time_rupture
    d['events'] = [events[0]]
    d['unbroken'] = [unbroken[0]]


def merge_extracted_values(extracted_merged, delta_values, rupture_values, events_values):
    """ save the extracted values to the merged datasets"""
    delta_ini, delta_fin = delta_values
    time_rupture, delta_rupture = rupture_values
    events, unbroken = events_values
    # ---
    deltas_merged, rupture_merged, events_merged = extracted_merged
    delta_ini_merged, delta_fin_merged = deltas_merged
    delta_rupture_merged, time_rupture_merged = rupture_merged
    molecules_with_events, unbroken_events = events_merged
    # --
    delta_ini_merged += delta_ini
    delta_fin_merged += delta_fin
    # --
    delta_rupture_merged += delta_rupture
    time_rupture_merged += time_rupture
    # --
    if events[0] > 0:
        molecules_with_events[0] += 1
    unbroken_events[0] += unbroken[0]
    return extracted_merged


def store_merged_extracted_values(d, extracted_merged):
    """ save the extracted values to the dictionary """
    deltas_merged, rupture_merged, events_merged = extracted_merged
    delta_ini_merged, delta_fin_merged = deltas_merged
    delta_rupture_merged, time_rupture_merged = rupture_merged
    molecules_with_events, unbroken_events = events_merged
    # --
    d['delta_ini_merged'] = delta_ini_merged
    d['delta_fin_merged'] = delta_fin_merged
    d['delta_rupture_merged'] = delta_rupture_merged
    d['time_rupture_merged'] = time_rupture_merged
    d['molecules_with_events'] = molecules_with_events[0]
    d['unbroken_events'] = unbroken_events[0]
