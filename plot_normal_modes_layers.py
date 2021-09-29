from matplotlib.pyplot import legend
import numpy as np
import sys
import os
import glob
from pathlib import Path
import matplotlib
from numpy.testing._private.utils import jiffies
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)
line_width = 0.4
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.6
matplotlib.rcParams['lines.linewidth'] = 0.6

markers = ['o','s','^','.','>','v']
colours = ['red','navy','mediumorchid','maroon','dodgerblue','gold']
linestyles = ['-','-','-','--','--','--']

x = []
y = []
results = {'geom' : [], 'idx' : [],'val' : [], 'co_freq' : [], 'sa_freq' : []}

normal_mode_outputs = sys.argv[1:]


for n,normal_mode_output in enumerate(normal_mode_outputs):
    parsing = False
    frequencies = []

    dir_name = os.path.dirname(normal_mode_output)
    n_layers = int(dir_name.split('/')[-1])

    print(normal_mode_output)

    with open(normal_mode_output,'r') as f:
        for line in f:
            if '-----' in line:
                continue
            if 'meV' in line:
                parsing = True
                continue
            if 'Zero-point energy' in line:
                parsing = False
                break

            if parsing:
                if 'i' in line:
                    without_imag = line.replace('i','')
                    frequencies.append(float(without_imag.split()[-1])) #cm-1
                else:
                    frequencies.append(float(line.split()[-1])) #cm-1

    data = np.array(frequencies)
    results['idx'].append(n_layers)
    results['val'].append(data)
    results['co_freq'].append(data[-1])
    results['sa_freq'].append(data[-2])



fig, ax = plt.subplots(1, 1)

idxs = np.array(results['idx'])

print(idxs)

#vals = np.array(results['val'])

#print(np.shape(vals))

x = idxs
y = np.array(results['co_freq'])

c=0
ax.plot(x,y,#label=r'$\Lambda_{{{}{}}}$'.format(labels[i],labels[j]),
    marker=markers[c],color=colours[c],linestyle=linestyles[c],
    mfc='none', markersize=4, label= 'IS')

c=1              
y = np.array(results['sa_freq'])
# ax.plot(x,y,#label=r'$\Lambda_{{{}{}}}$'.format(labels[i],labels[j]),
#     marker=markers[c],color=colours[c],linestyle=linestyles[c],
#     mfc='none', markersize=4, label = 'SA')


# ax.set_ylabel(r'$\Lambda_{\mathrm{ij}}$ / ps$^{-1}$',color='black')
ax.set_ylabel(r'$\omega$ / cm$^{-1}$',color='black')
ax.set_xlabel(r'$N_{\mathrm{layers}}$',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

ax.legend(fancybox=True,framealpha=0)
#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.15), loc='center')
fig.set_figheight(4)
fig.set_figwidth(5)
fig.savefig('normal_mode_layers.pdf',transparent=True,bbox_inches='tight')
