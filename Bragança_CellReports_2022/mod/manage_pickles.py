import os
import pickle


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


def save_data_to_pickle(info, area, data_to_store):
    _, condition, date, flowcell = info
    new_exp = f'{date}_cell#{flowcell}_{area}'
    # ---
    filepath = os.path.join('pickles', f'{condition}.pickle')
    if os.path.exists(filepath):
        # Load data, deserialize saved dictionary
        print(f'Loading stored dictionary...', end='')
        with open(filepath, 'rb') as file_handler:
            stored = pickle.load(file_handler)
        print(list(stored.keys()))
        # Add new data to the existing dictionary
        stored[new_exp] = data_to_store
        data_to_file = stored
    else:
        print('Creating new dictionary.')
        data_to_file = {new_exp: data_to_store}
        
    # save to pickle
    with open(filepath, 'wb') as file_handler:
        pickle.dump(data_to_file, file_handler)
    file_handler.close()
    print(f'Data from {new_exp} saved.')