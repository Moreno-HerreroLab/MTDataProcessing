import os
import pandas as pd


def import_experimental_details(folder):
    path = os.path.join(folder,'experimental_details.txt')
    if os.path.exists(path): details = pd.read_csv(path, sep='\t')
    else: details = pd.DataFrame(columns=['condition', 'date', 'flowcell#', 'area#', 'results#', 'total_sample', 'discard#', 'final_sample'])
    return details


def save_experimental_details(info, data, line, details):
    folder, condition, date, flowcell = info
    area, R, del_beads = line
    total = data[area]['total_sample']
    valid = data[area]['valid_sample, N']
    new_row = { 'condition': condition, 'date': date, 'flowcell#': [flowcell],
                'area#': [area], 'results#': [R], 'total_sample': [total], 
                'discard#': [del_beads], 'final_sample': [valid]}
    # --- is there any previous stoeed info for that experiment?
    cond1 = (details['date'] == date)
    cond2 = (details['flowcell#'] == flowcell)
    cond3 = (details['area#'] == area)
    idx = details.index[cond1 & cond2 & cond3]
    # --- no previous info, create new row
    if idx.empty:
        new_row = pd.DataFrame(new_row)
        details = pd.concat([details, new_row], axis=0, ignore_index=True)
    # --- yes, previous stored info, replace it
    else: 
        details.loc[idx[0]:idx[0], :] = new_row.values()
    # --- save
    path = os.path.join(folder,'experimental_details.txt')
    details.to_csv(path, index=None, sep='\t')
    return details
    