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
tensors = sys.argv[1:]
x = []
y = []
results = {'geom' : [], 'idx' : [],'val' : []}


labels = ['\mathrm{C}_x','\mathrm{C}_y','\mathrm{C}_z','\mathrm{O}_x','\mathrm{O}_y','\mathrm{O}_z']

for t,tensor in enumerate(tensors):
    print(tensor)
    data = np.loadtxt(tensor)
    dir_name = os.path.dirname(tensor).replace('../../../opt/','')
    n_layers = int(dir_name.split('/')[-1])

    if n_layers > 16:
        flipped_data = np.zeros_like((data))
        flipped_data[0:3,0:3]=data[3:,3:]
        flipped_data[3:,3:]=data[0:3,0:3]
        results['val'].append(flipped_data)
    else:
        results['val'].append(data)

    results['idx'].append(n_layers)



fig, ax = plt.subplots(1, 1)

idxs = np.array(results['idx'])

print(idxs)

vals = np.array(results['val'])
dimension = np.shape(vals)[1]
c=0
for i in range(dimension):
    #if i in [0,5]:
        for j in range(dimension):
            if i == j: #or j == 3:
                x = idxs
                y = vals[:,i,j]
                ax.plot(x,y,label=r'$\Lambda_{{{}{}}}$'.format(labels[i],labels[j]),marker=markers[c],color=colours[c],linestyle=linestyles[c],
                mfc='none', markersize=4)
                c += 1


ax.set_ylabel(r'$\Lambda_{\mathrm{ij}}$ / ps$^{-1}$',color='black')
ax.set_xlabel(r'$N_{\mathrm{layers}}$',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))

ax.set_ylim(bottom=0.0)
ax.legend(fancybox=True,framealpha=0)
#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.15), loc='center')
fig.set_figheight(4)
fig.set_figwidth(5)
fig.savefig('fig.pdf',transparent=True,bbox_inches='tight')
