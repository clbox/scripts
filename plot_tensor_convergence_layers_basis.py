import numpy as np
import sys
import os
import glob
from pathlib import Path
import matplotlib
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
    results['val'].append(data)
    results['basis'].append(basis)



fig, ax = plt.subplots(1, 1)

idxs = np.array(results['idx'])
vals = np.array(results['val'])
dimension = np.shape(vals)[1]

list_of_used_basis =  np.array(results['basis'])
for b,basis in enumerate(['light+','tight','really_tight']):

    basis_index = np.argwhere(list_of_used_basis==basis)
    c=0
    for i in range(dimension):
        #if i in [0,5]:
            for j in range(dimension):
                #if i == j: #or j == 3:
                if i == 2 and j == 2:
                    x = idxs[basis_index]
                    y = vals[basis_index,i,j]
                    # ax.plot(x,y,label=r'$\Lambda_{{{}{}}}$'.format(labels[i],labels[j]),marker=markers[c],color=colours[c],linestyle=linestyles[c],
                    ax.plot(x,y,label=basis,marker=markers[b],color=colours[b],linestyle=linestyles[b],
                    mfc='none', markersize=4)
                    c += 1


ax.set_ylabel(r'$\Lambda_{\mathrm{ij}}$ / ps$^{-1}$',color='black')
ax.set_xlabel(r'$N_{\mathrm{layers}}$',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))

ax.legend(fancybox=True,framealpha=0)
#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.15), loc='center')
fig.set_figheight(4)
fig.set_figwidth(5)
fig.savefig('fig.pdf',transparent=True,bbox_inches='tight')
