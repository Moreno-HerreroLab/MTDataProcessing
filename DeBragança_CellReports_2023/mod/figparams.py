import matplotlib.pyplot as plt

def def_fig_parameters():
    cm = 1/2.54 # conversion: 1 centimeter in inches
    plt.rc('figure', titlesize=12, figsize=(15*cm, 10*cm))
    plt.rc('axes', titlesize=12, labelsize=11)  # font size of the axes title and of the x and y labels
    plt.rc('xtick', labelsize=10, direction='in', top=True)  # font size of the tick labels
    plt.rc('ytick', labelsize=10, direction='in', right=True)  # font size of the tick labels
    plt.rc('legend', fontsize=10, markerscale=10, loc='lower left')  # legend font size
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    print('Figure parameters adjusted.')
    return cm