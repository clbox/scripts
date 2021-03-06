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

style = sys.argv[1] #k (k grid), b (broadening), rc (reaction coordinate), cov (coverage), unit (unit cell), o (other)
markers = ['o','s','^','.','>','v']
colours = ['red','navy','mediumorchid','maroon','dodgerblue','gold']
linestyles = ['-','-.','--','-','-.','--']

tensors = sys.argv[2:]
x = []
y = []
results = {'geom' : [], 'idx' : [],'val' : []}

for t,tensor in enumerate(tensors):
    print(tensor)
    data = np.loadtxt(tensor)



    results['idx'].append(os.path.dirname(tensor))
    results['val'].append(data)



fig, ax = plt.subplots(1, 1)

idxs = np.array(results['idx'])
vals = np.array(results['val'])
dimension = np.shape(vals)[1]
c=0
for i in range(dimension):
    for j in range(dimension):
        if i == j:
            x = idxs
            if style == 'rc':
                x = np.arange(0,(len(idxs)))
            if style == 'unit':
                x = [r'$2 \times 2 \times 4$',r'$3 \times 3 \times 4$',r'$4 \times 4 \times 4$']
            y = vals[:,i,j]
            ax.plot(x,y,label=r'$\Lambda_{{{}{}}}$'.format(i+1,j+1),marker=markers[c],color=colours[c],linestyle=linestyles[c],
            markersize=3,linewidth=1)
            c += 1


ax.set_ylabel(r'$\Lambda_{\mathrm{ij}}$ / ps$^{-1}$',color='black')
if style =='k':
    ax.set_xlabel(r'$N_{\mathrm{k}}$',color='black')
elif style =='rc':
    ax.set_xlabel(r'Reaction coordinate',color='black')
elif style =='unit':
    ax.set_xlabel(r'Unit cell size',color='black')
else:
    ax.set_xlabel(r'Index',color='black')


ytop = np.max(vals)
ytop = np.ceil(ytop)

ymin = np.min(np.diagonal(vals,axis1=1,axis2=2))
ymin = np.floor(ymin)

ax.set_ylim(ymin,ytop)

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.2))

if style == 'rc':
    if len(tensors) > 10:
        ax.xaxis.set_major_locator(MultipleLocator(2))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(1))

    if ytop > 5 and ymin == 0:
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(1.0))

    



plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.15), loc='center')
fig.set_figheight(2.)
fig.set_figwidth(3.25)
fig.savefig('fig.pdf',transparent=True,bbox_inches='tight')