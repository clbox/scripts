from matplotlib.pyplot import legend
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


markers = ['o','s','^']
colours = ['navy','forestgreen','maroon']

legend_labels = ['Adsorption', 'Reactant','Transition state']

i = int(sys.argv[1])
j = int(sys.argv[2])
tensors = sys.argv[3:]
x = []
y = []
results = {'geom' : [], 'cov' : [],'val' : []}

for t,tensor in enumerate(tensors):
    data = np.loadtxt(tensor)


    results['geom'].append(os.path.dirname(os.path.dirname(tensor)))
    results['cov'].append(int(os.path.dirname(tensor)[-1]))
    results['val'].append(data[i,j])



fig, ax = plt.subplots(1, 1)

geoms = np.array(results['geom'])
covs = np.array(results['cov'])
vals = np.array(results['val'])

unique = list(set(results['geom']))
unique.sort()
for g,geom in enumerate(unique):
    print(geom)
    idx = np.argwhere((geoms==geom)).flatten()
    x = covs[idx]
    y = vals[idx]
    ax.plot(x,y,label=geom,marker=markers[g],linestyle='-',color=colours[g])


ax.set_ylabel(r'$\Lambda_{\mathrm{rr}}$ / ps$^{-1}$',color='black')
ax.set_xlabel(r'Unit cell',color='black')

ticks=[]
tick_labels = []
labels= set(results['cov'])
for i in labels:
    ticks.append(int(i))
    tick_labels.append(r'p({}$\times${})'.format(i,i))

ax.set_xticks(ticks)
ax.set_xticklabels(tick_labels)
ax.set_ylim(0,1)



handles,labels = ax.get_legend_handles_labels()

handles = handles
labels = legend_labels

ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.1))


plt.legend(handles=handles,labels=labels,ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.1), loc='center')
fig.set_figheight(2.)
fig.set_figwidth(3.25)
fig.savefig('fig.pdf',transparent=True,bbox_inches='tight')