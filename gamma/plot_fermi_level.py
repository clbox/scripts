import gamma.friction_coupling as fc
import glob
import sys
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

line_width = 0.8
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.8
matplotlib.rcParams['lines.linewidth'] = 1.2

fig, ax = plt.subplots(1, 1)





dir_names = sys.argv[1:]

dir_names.sort()

print(dir_names)

time_step = 1 #fs
time_axis = np.arange(0,len(dir_names))
time_axis = time_axis * time_step / 1000


chem_pots = []
parser = fc.friction_output_parser_2021()
for dir_name in dir_names:
    gamma_files = glob.glob(dir_name+'/*friction_gamma*.out')
    chem_pot = parser.parse_chem_pot(dir_name+'/aims.out')
    print('Chemical potential / eV : '+ str(chem_pot))
    chem_pots.append(chem_pot)



ax.plot(time_axis,chem_pots,color='black',linestyle='none',mfc='orangered',marker='o')



ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.5))

# ax.yaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(1))

ax.set_xlim(0.,2.4)

ax.set_xlabel(r'$t$ / $\mathrm{ps}$')
ax.set_ylabel(r'$\epsilon_\mathrm{F}$ / $\mathrm{eV}$')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')



fig.set_figheight(2)
fig.set_figwidth(3)

fig.savefig('fermi_level.pdf',transparent=False,bbox_inches='tight')