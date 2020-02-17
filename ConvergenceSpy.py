import sys, os, subprocess

import numpy as np
from numpy.ma import masked_where as mask

import matplotlib.pyplot as plt
from matplotlib import ticker

from ase.io import read

"""
This script will track the convergence proprieties of a given FHI-aims
calculation as a function of SCF-cycle number.
The script is invoked using the following command:
    - python ConvergenceSpy.py /path/to/aims/calc/dir
With the optional arguments:
    - log: This will use logarithmic axes.
    - save: This will save the image to a .png file.

Converged & unconverged points are marked with green & red dots respectively.
Settings which have converged have their plots outlined in green, those that
have not are outlined in red. The convergence window is marked by two blue
horizontal lines, points which land between them are said to be converged.
"""
def show_plot(fhi_aims_dir):

    def set_plot(path):


        # Read in the conv.dat file and get the last SCF cycle
        lines = open(f'{path}/conv.dat', 'r').readlines()
        last_SCF = lines[-int(lines[-1].split()[0]):]

        # Retrieve last SCF cycle's data. The ignored value is time taken for
        # each step.[SCF, D_t, D_rho, D_EV (eV), D_E (eV), D_F (eV/A)]
        step, _, rho, EV, E, F = np.array([l.split()[:6] for l in last_SCF]).T

        # This lists which convergence criteria were met at each step
        conv_step = np.asarray([['A_rho' in l, 'A_EV' in l, 'A_E' in l, 'A_F' in l]
                     for l in last_SCF]).T

        # Fetch the convergence settings for rho, EV, E, and F (in order)
        cov_set = list(map(float, lines[5][1:].split()))

        #___Data Handling___#
        # Convert the data into floats and integers as needed
        rho, EV, E = rho.astype('float'), EV.astype('float') , E.astype('float')
        step = step.astype('int')

        try:  # F is handled differently as its "[not calc.]" in the 1st SCF
            F = F.astype('float')
        except ValueError:  # If no force, state it on the plot
            F = np.ones(len(step)) * 1E-4
            axs.flat[3].text(0.3, 0.5, 'Not Calculated', va='center',
                             transform=axs.flat[3].transAxes,
                             bbox={'facecolor': 'red'})

        #___Plot Setup___#
        # If using logarithmic scale, take the absolute of all values
        if LOG_SCALE:
            rho, EV, E, F = abs(np.array([rho, EV, E, F]))

        plot_data = zip(
            # Subplots
            axs.flat,

            # Y-axis data to be plotted on said subplots
            [rho, EV, E, F],

            # Titles of the subplots
            [rf'|$\Delta$ {i}|' for i in ['Density', 'Eigen Energy',
                                          'Total Energy', 'Forces']],

            # Labels for the subplots' axes
            [[None, r'$\rho$'], [None, r'E$_{Eigen}$ (eV)'],
             ['Step', r'E$_{Total}$ (eV)'], ['Step', r'Force (eV/$\AA$)']],

            conv_step, cov_set)

        #___Plot Building___#
        for ax, ys, title, label, converged, thresh in plot_data:
            # Black line with red/green dots on unconverged/converged points.
            ax.plot(step, ys, '-k',
                    step, mask(converged, ys), '.r', 
                    step, mask(np.invert(converged), ys), '.g')

            ax.axhline(y=thresh) # Convergence point line
            ax.set_title(title)  # Subplot title
            ax.set(xlabel=label[0], ylabel=label[1])  # Set subplot axes

            #___Scale Setup___#
            # If using logarithmic scale, set the scale to Log10
            if LOG_SCALE:
                ax.set_yscale('log')
            else:  # If not, then:
                # Use scientific notation
                fmt = ticker.ScalarFormatter(useMathText=True) 
                fmt.set_scientific(True) 
                fmt.set_powerlimits((-1,1))
                ax.yaxis.set_major_formatter(fmt)

                # Rescale the y-axis ignoring the 1st and last values,
                # as then tend to be too big and destroy the scale.
                y_max, y_min = max(ys[1:-1]), min(ys[1:-1])
                stretch = abs(y_max - y_min)/10
                ax.set_ylim([y_min - stretch, y_max + stretch])

                # As we have positive and negative numbers, we need to add a
                # second convergence threshold line on the negative side
                ax.axhline(y=-thresh)

            # Colour the plot's box green or red based on convergence
            colour = {True: 'green', False: 'red'}[converged[-1]]
            for i in ['bottom','top','right','left']:
                ax.spines[i].set_color(colour)
                ax.spines[i].set_linewidth(1.2)

            # Add in the threshold value
            ax.text(0.98, 0.98, f'Threshold: {thresh:7.1e}', va='top',
                    ha='right', transform=ax.transAxes, fontdict={'size': 8})

        # Get the chemical formula of the system and set the plot's title
        system = read(f'{path}/geometry.in').get_chemical_formula()
        fig.suptitle(f'Convergence plots for system {system} ({os.path.basename(path)})')

    # Construct a 2x2 plot
    fig, axs = plt.subplots(2, 2, sharex=True, constrained_layout=True)
    set_plot(fhi_aims_dir)
    plt.show()

    if SAVEFIG:
        plt.savefig(f'{os.path.basename(fhi_aims_dir)}.png',
                    bbox_inches='tight')


global LOG_SCALE
if 'log' in sys.argv:
    LOG_SCALE = True
else:
    LOG_SCALE = False


global SAVEFIG
if 'save' in sys.argv:
    SAVEFIG = True
else:
    SAVEFIG = False

aims_dir=sys.argv[1]
show_plot(aims_dir)
