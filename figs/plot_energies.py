import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
import glob
from mem_energy_loss import read_memory_kernel

SMALL_SIZE = 9.5
MEDIUM_SIZE = 9.5
BIGGER_SIZE = 9.5

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
annotate_args = {'xy' : (0.05,0.9), 'xycoords' : 'axes fraction'}

filename = sys.argv[1]

pathway = np.loadtxt(filename)

min_e = np.min(pathway[:,1])

pathway[:,1] = pathway[:,1] - min_e

barrier = np.max(pathway[:,1]) - np.min(pathway[:,1])

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')



ax.plot(pathway[:,0],pathway[:,1],marker='s',mfc='None',linestyle='--',color='grey',markersize=4,linewidth=1)
ax.annotate(r'$\Delta E = {:.2f}$ eV'.format(barrier),ha="left", **annotate_args)


ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(0.025))
ax.yaxis.set_major_locator(MultipleLocator(0.05))

ax.set_ylim(bottom=0)

ax.set_xlabel("Reaction coordinate")
ax.set_ylabel(r'Energy / eV')
fig.set_figheight(2.0)
fig.set_figwidth(3.25)
fig.savefig('fig.pdf',transparent=True,bbox_inches='tight')