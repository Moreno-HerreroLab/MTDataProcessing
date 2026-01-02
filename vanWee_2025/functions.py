import matplotlib.pyplot as plt
import numpy as np
import os
import gc

# ----------------------------------------

# Relation between the distance to the surface ("shift pos", or magnet) in the data and the applied force.
# This equation depends on the magnet configuration. 
force_equation = lambda magnet: 5.66*np.exp(-0.68*magnet)   #(2025_Mag5_singlelayer_T1-2024)


def add_results_to_experiments_log(experiments_log, condition, results):
    """
    Add results from one condition to the global experiment log.
    """
    for line in results:
        line.insert(0, condition)
        experiments_log.append(line)

    return experiments_log


def process_step(z_defined_step, time_defined_step, 
                 diff_threshold_um, frames_threshold, rolling_window):
    """
    Process individual step:
    - smooth the data
    - detect plateaux
    - compute lifetimes
    """
    # Smooth the signal with a moving average
    kernel = np.ones(rolling_window) / rolling_window
    smoothed = np.convolve(z_defined_step, kernel, mode="same")

    # Find plateaux in the smoothed signal
    plateaux, diff, zplateaux = find_plateaux(
        smoothed, diff_threshold_um, frames_threshold)

    # Compute lifetimes of the plateaux
    lifetimes = extract_lifetimes(time_defined_step, plateaux)
    filtered = [(z, 1) if z >= -0.05 else (z, t) for (z, t) in zip(zplateaux, lifetimes)]
    zplateaux, lifetimes = zip(*filtered) if filtered else ([], [])
    zplateaux, lifetimes = list(zplateaux), list(lifetimes)

    return smoothed, diff, zplateaux, lifetimes


def find_plateaux(signal, diff_threshold_um, frames_threshold):
    """
    Detect flat regions (plateaux) in a signal.
    """
    # Compute the absolute difference between consecutive elements
    diff = np.abs(np.diff(signal))
    # Identify continuous plateaux
    flat_regions = diff < diff_threshold_um
    plateaux = []
    zplateaux = []
    start = None
    for i, is_flat in enumerate(flat_regions):
        # is_flat == True ---> plateaux
        if is_flat and start is None: start = i  # start a new plateau
        # is_flat == False ---> rupture!!
        elif not is_flat and start is not None: # Evaluate the detected plateaux 
            if (i-start) >= frames_threshold: # if it's long enough, save it
                plateaux.append((start, i))
                delta = np.mean(signal[start:i])
                zplateaux.append(delta)
            else: # ignore it if it's too small 
                delta = None
                pass  
            start = None # restart
            if delta != None and delta >= -0.05: break
    # Add the last plateau
    if start is not None and (i-start) >= frames_threshold:
        plateaux.append((start, i))
        zplateaux.append(np.mean(signal[start:i]))
    return plateaux, diff, zplateaux



def extract_lifetimes(time_defined_step, plateaux):
    """
    Compute normalized lifetimes for each plateau.
    """
    lifetimes = []
    ini = time_defined_step.iloc[0]
    fin = time_defined_step.iloc[-1]
    for start, end in plateaux:
        lifetime = (time_defined_step.iloc[end] - ini) / (fin - ini)
        lifetimes.append(lifetime)

    return lifetimes


def plot_and_export_individual_events_with_details(z, time, force, steps, label, diff_threshold_um, time_minimum_ms, rolling_window, figpath):
    frame_rate = 120  # experimental value, frames/s
    frames_threshold = int(frame_rate * (time_minimum_ms*10**(-3)))

    fig, ax = plt.subplots(figsize=(6,4), nrows=3, ncols=2, sharey='row', sharex='col',
                            gridspec_kw={'width_ratios': [5, 1], 'height_ratios': [4, 1, 1],
                                        'top': .91, 'bottom': .13, 'hspace': 0.025,
                                        'wspace': 0.0115, 'left': .15, 'right': .99})

    # define list of axis of interest and remove unused axis:
    axlist = (ax[0,0], ax[1,0], ax[0,1], ax[2,0])
    for axis in (ax[1,1], ax[2,1]): axis.remove()

    for step in range(2, len(steps.unique())+1, 2):
        # step data
        time_defined_step = time[steps==step]
        z_defined_step = z[steps==step]
        smoothed, diff, zplateaux, lifetimes = process_step(z_defined_step, time_defined_step, diff_threshold_um, frames_threshold, rolling_window)

        # Raw and smoothed data
        ax[0,0].scatter(time_defined_step, z_defined_step, s=0.5, color='lightgrey', alpha=1, label='raw data')
        ax[0,0].plot(time_defined_step, smoothed, color='black', alpha=1, label=f'Convolution (kernel window)')
        ax[0,0].set(ylabel='DNA reduction, Δz (µm)', ylim=(-1.59, 0.19))
        
        # Plot plateaux
        start = time_defined_step.min()
        end = time_defined_step.max()
        for n, (z_plateau, life) in enumerate(zip(zplateaux, lifetimes)):
            t_plateau = (end-start)*life + start
            ax[0,0].plot([start, t_plateau], [z_plateau, z_plateau], color='red', lw=2, ls=':', 
                         label=f'Plateaus (threshold={diff_threshold_um}nm, time_limit={time_minimum_ms}ms ({frames_threshold}frames)' if n==0 else None)
            ax[0,0].annotate(text=f'Δz = {round(z_plateau, 2)}', xy=(t_plateau-0.7, z_plateau-0.10), va='top', color='red', rotation=90, fontsize=8)
            ax[0,0].annotate(text=f'tau = {round(life, 2)}', xy=(t_plateau-0.4, z_plateau-0.10), va='top', color='red', rotation=90, fontsize=8)
        ax[0,0].legend(loc='lower right', fontsize=8)

        # Plot diff
        ax[1,0].plot(time_defined_step[:-1], diff*1000, color='grey', lw=1)
        ax[1,0].axhline(y=diff_threshold_um, color='red', linestyle=':', linewidth=1)
        ax[1,0].set(ylabel='Diff (nm)', ylim=(0, 15))
        
        # Plot force
        ax[2,0].plot(time_defined_step, force[steps==step], color='grey', lw=1)
        ax[2,0].set( xlabel='Time (s)', xlim=(start, end), ylabel='F (pN)', ylim=(-0.1, 3.1))
        
        # Plot histogram
        bins = np.arange(min(z_defined_step), max(z_defined_step), 0.03)
        ax[0,1].hist(z_defined_step, bins=bins, color='grey', edgecolor='white', linewidth=1, orientation='horizontal')
        ax[0,1].set(xlabel='Density', xlim=(0, None))
        
        # Grid styling
        for axis in axlist:
            axis.grid(color='lightgrey', linestyle=':', linewidth=1)
            axis.set_axisbelow(True)

        # Update plots and export
        fig.savefig(os.path.join(figpath, f'{label}_Step#{step}.png'), transparent=False, dpi=150)

        # Clean up memory for each iteration
        del smoothed
        for axis in axlist: axis.clear()
        gc.collect()  # Force garbage collection to free memory
        plt.close(fig)

