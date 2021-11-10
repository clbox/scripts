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
matplotlib.rcParams['lines.markeredgewidth'] = 0.5
matplotlib.rcParams['lines.linewidth'] = 0.8


markers = ['s','o','^','.','>','v']
colours = ['orangered','dodgerblue','seagreen','mediumorchid','maroon','dodgerblue','gold']
linestyles = ['-','-','-','--','--','--']
tensors = sys.argv[1:]
x = []
y = []
results = {'geom' : [], 'idx' : [],'val' : [], 'basis' : []}


labels = ['\mathrm{C}_x','\mathrm{C}_y','\mathrm{C}_z','\mathrm{O}_x','\mathrm{O}_y','\mathrm{O}_z']

for t,tensor in enumerate(tensors):
    print(tensor)
    data = np.loadtxt(tensor)
    dir_name = os.path.dirname(tensor)
    basis = dir_name.split('/')[0]
    n_layers = int(dir_name.split('/')[-1])
    results['idx'].append(n_layers)

    if 'second_batch' in tensor:
        data_adjust = np.zeros((6,6))
        data_adjust[0:3,0:3] = data[3:6,3:6]
        results['val'].append(data_adjust)
    else:
        results['val'].append(data)
    results['basis'].append(basis)



fig, ax = plt.subplots(1, 1)

idxs = np.array(results['idx'])
vals = np.array(results['val'])
dimension = np.shape(vals)[1]

list_of_used_basis =  np.array(results['basis'])
for b,basis in enumerate(['hgt','c2']):

    basis_index = np.argwhere(list_of_used_basis==basis)
    c=0
    for i in range(dimension):
        #if i in [0,5]:
            for j in range(dimension):
                #if i == j: #or j == 3:
                if i == 2 and j == 2:
                    if basis == 'c2':
                        label = r'$\boldsymbol{G}$'
                    else:
                        label= r'$\boldsymbol{G}^{\mathrm{HGT}}$'


                    x = idxs[basis_index]
                    y = vals[basis_index,i,j]
                    # ax.plot(x,y,label=r'$\Lambda_{{{}{}}}$'.format(labels[i],labels[j]),marker=markers[c],color=colours[c],linestyle=linestyles[c],
                    ax.plot(x,y,label=label,marker=markers[b],color=colours[b],linestyle='none',
                    mfc=colours[b], markeredgecolor='black', markersize=4)
                    c += 1


#ax.set_ylabel(r'$\Lambda_{\mathrm{ij}}$ / ps$^{-1}$',color='black')
ax.set_ylabel(r'$\Lambda_{\mathrm{C}_z,\mathrm{C}_z}$ / meV ps$^{-1}$ $\mathrm{\AA{}}^{-2}$',color='black')
ax.set_xlabel(r'$N_{\mathrm{layers}}$',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.xaxis.set_major_locator(MultipleLocator(4))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))
ax.yaxis.set_major_locator(MultipleLocator(0.05))

ax.set_xlim(3,26)
ax.set_ylim(0.1,0.35)

ax.legend(fancybox=True,framealpha=1,edgecolor='black',handletextpad=0.05,borderpad=0.3,handlelength=1.2)

#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.15), loc='center')
fig.set_figheight(3)
fig.set_figwidth(2)
fig.savefig('coupling.pdf',transparent=True,bbox_inches='tight')
