import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import shutil

def def_fig_parameters():
    plt.rc('axes', titlesize=14, labelsize=12)  # font size of the axes title and of the x and y labels
    plt.rc('xtick', labelsize=10, direction='in', top=True)  # font size of the tick labels
    plt.rc('ytick', labelsize=10, direction='in', right=True)  # font size of the tick labels
    plt.rc('legend', fontsize=10, markerscale=10, loc='lower left')  # legend font size
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

def check_directories(input_folder):
    corr_folder = os.path.join(input_folder, 'data_corrected')
    if not os.path.exists(corr_folder): os.mkdir(corr_folder)
    ori_folder = os.path.join(input_folder, 'data_original')
    if not os.path.exists(ori_folder): os.mkdir(ori_folder)
    figs_folder = os.path.join(corr_folder, 'allBeadsCorrected')
    if not os.path.exists(figs_folder): os.mkdir(figs_folder)
    return ori_folder, corr_folder, figs_folder


def import_file(input_folder, file):
    PAR = False
    # -- import --
    filename = f'resultados_{file}.dat'
    filepath = os.path.join(input_folder, filename)
    if os.path.exists(filepath):
        data = pd.read_csv(filepath, delimiter='\t', decimal=',', header=0)
        data = data.drop(data.columns[-1], axis=1)
    else:
        filename = f'resultados_{file}_PAR.dat'
        filepath = os.path.join(input_folder, filename)
        error_message = f'No resultados_{file}.dat nor resultados_{file}_PAR.dat file was found.'
        assert (os.path.exists(filepath)), f'{error_message}'
        PAR = True
        data = pd.read_csv(filepath, delimiter='\t', header=0)
    return filename, data, PAR


def get_labels(data):
    t_label = [col for col in data.columns if 'Time' in col]
    xyz_labels = [col for col in data.columns if 'Bead' in col]
    x_labels = [col for col in xyz_labels if 'X' in col]
    y_labels = [col for col in xyz_labels if 'Y' in col]
    z_labels = [col for col in xyz_labels if 'Z' in col]
    num = int(len(z_labels))
    names = [label.split('-')[0] for label in z_labels]
    return names, num, x_labels, y_labels, z_labels, t_label


def check_beads_num(data, num, z0_file, exp_file):
    tethers = [col for col in data.columns if 'Bead' in col and 'Z' in col]
    tethers = len(tethers)
    if num != tethers:
        print(f'{z0_file} has {num} beads and {exp_file} has {tethers} beads.')
        print('Some beads will raise an error and will be ignored.')


def move_original_data(filename, input_folder, ori_folder):
    current_path = os.path.join(input_folder, filename)
    destination_path = os.path.join(ori_folder, filename)
    shutil.move(current_path, destination_path)


def z0correction(input_folder, z0_file, id_file, exp_file, showRaw=False, DNA_length=2.5, yaxlim=None):
    yaxlim = (-0.5, DNA_length + 1) if yaxlim is None else yaxlim
    ori_folder, corr_folder, figs_folder = check_directories(input_folder)
    # ----- import z0 -----
    filename, data_z0, _ = import_file(input_folder, z0_file)
    move_original_data(filename, input_folder, ori_folder)
    # --------- labels ------
    names, num, x_labels, y_labels, z_labels, t_label = get_labels(data_z0)
    # --------- z_min --------
    zmin = [data_z0[z_labels[i]].min() for i in range(0, num)]
    # -- Plot z0_file for visual validation
    for i in range(num):
        t = data_z0[t_label[0]]
        x = data_z0[x_labels[i]]
        y = data_z0[y_labels[i]]
        z = data_z0[z_labels[i]]

        fig, axs = plt.subplots(nrows=3, sharex=True, gridspec_kw={'height_ratios': [3, 1, 1]})
        fig.suptitle(f'results_{z0_file}_z0_' + names[i])
        axs[0].set_title('z position')
        axs[0].scatter(t, z, marker='o', color='blue', s=1)
        axs[0].set(ylabel='z (um)', ylim=yaxlim)
        axs[0].annotate(f'zmin = {zmin[i]} um', xy=(1, -0.5), color='red')
        axs[0].axhline(zmin[i], 0, 1, ls=':', color='red', ms=1)

        axs[1].set_title('x position')
        axs[1].scatter(t, x, marker='o', color='blue', s=1)
        axs[1].set(ylabel='x (um)', ylim=(-1, 1))

        axs[2].set_title('y position')
        axs[2].scatter(t, y, marker='o', color='blue', s=1)
        axs[2].set(xlabel='Time (s)', ylabel='y (um)')
        axs[2].set(xlim=(0, t.max()), ylim=(-1, 1))

        figname = f'results_{z0_file}_z0_' + names[i] + '.png'
        figpath = os.path.join(figs_folder, figname)
        plt.savefig(figpath, facecolor='w')
        plt.close(fig)

    # ------ function to correct files ----
    def correction(i, zmin, data, data_corr, t_label, names, file, figs_folder, DNA_length, ID, PAR, showRaw):
        # -- substract zmin --
        data_corr[z_labels[i]] = data_corr[z_labels[i]] - zmin[i]
        # ---- plot for visual validation ----
        t = data[t_label[0]]
        x = data[x_labels[i]]
        y = data[y_labels[i]]
        z = data[z_labels[i]]
        zcorr = data_corr[z_labels[i]]

        fig, axs = plt.subplots(nrows=3, sharex=True, gridspec_kw={'height_ratios': [3, 1, 1]})
        fig.suptitle(f'results_{file}_' + names[i])

        axs[0].set(title='z position')
        axs[0].scatter(t, zcorr, marker='o', color='orange', s=1, label='corr data')
        axs[0].set(ylabel='z (um)', xlim=(0, x.max()), ylim=yaxlim)
        if showRaw:
            axs[0].scatter(x, z, marker='o', color='green', s=1, label='raw data')
            axs[0].legend(markerscale=10)
        if ID:
            figname = f'results_{file}_ID_{names[i]}.png'
            full_Z = np.mean(zcorr[0:121])
            axs[0].axhline(full_Z, 0, 1, ls=':', color='black', ms=1)

        axs[1].set(title='x position')
        axs[1].scatter(t, x, marker='o', color='orange', s=1)
        axs[1].set(ylabel='x (um)', ylim=(-1, 1))
        axs[1].axhline(0.0, 0, 1, ls=':', color='black', ms=1)

        axs[2].set(title='y position')
        axs[2].scatter(t, y, marker='o', color='orange', s=1)
        axs[2].axhline(0.0, 0, 1, ls=':', color='black', ms=1)
        axs[2].set(xlabel='Time (s)', ylabel='y (um)', xlim=(0, t.max()), ylim=(-1, 1))

        if PAR:
            figname = f'results_{file}_PAR_{names[i]}.png'
        elif not ID:
            figname = f'results_{file}_{names[i]}.png'
        figpath = os.path.join(figs_folder, figname)
        plt.savefig(figpath, facecolor='w')
        plt.close(fig)
        return data_corr

    # ----- correct tethers-ID file ----
    ID = True
    PAR = False
    filename, data, _ = import_file(input_folder, id_file)
    data_corr = data.copy(deep=True)
    for i in range(num): correction(i, zmin, data, data_corr, t_label, names, id_file, figs_folder, DNA_length, ID, PAR, showRaw)
    output_file = f'results_{id_file}_ID_corr.dat'
    filepath = os.path.join(corr_folder, output_file)
    data_corr.to_csv(filepath, index=False, sep='\t')
    move_original_data(filename, input_folder, ori_folder)
    # ----- correct exp_files -----
    ID = False
    if type(exp_file) is int: exp_file = [exp_file]
    for file in exp_file:
        filename, data, PAR = import_file(input_folder, file)
        data_corr = data.copy(deep=True)
        # -- evaluate if z0 and exp-files have the same number of beads --
        check_beads_num(data, num, z0_file, file)
        # -- treat data & plot for validation --
        for i in range(num):
            try:
                correction(i, zmin, data, data_corr, t_label, names, file, figs_folder, DNA_length, ID, PAR, showRaw)
            except:
                print(f'Error! Problem with {names[i]} from {filename}.')
        if PAR:
            output_file = f'results_{file}_PAR_corr.dat'
        else:
            output_file = f'results_{file}_corr.dat'
        filepath = os.path.join(corr_folder, output_file)
        data_corr.to_csv(filepath, index=False, sep='\t')
        move_original_data(filename, input_folder, ori_folder)
    # ----------------------
    print('Done.')
    

def correction_by_areas(folder, z0, ID, R, DNA_length=1.5):
    assert len(z0)==len(ID) & len(ID)==len(R), 'The arrays are not the same size.'
    for i in range(0, len(z0)): 
        print(f'area {i}...', end=' ')
        z0correction(input_folder = folder, z0_file = z0[i], id_file = ID[i],
                          exp_file = R[i], DNA_length=DNA_length)
    print('Finished!')   