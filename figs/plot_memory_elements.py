#Ideally takes element number as first argument then a list of filenames of *_memory_kernel.out 

#Additionally can read in tensors and mark where the tensor evaluated to.

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
#labels = ['Reactant','Adsorption','Transition state','Dissociated']
#labels = ['Adsorption','Transition state','Dissociated']

#labels = list(range(0,15))

colours = ['lightgrey','grey','black']

labels =[]
# for i in [2,3,4]:
#     labels.append(r'p({}$\times${})'.format(i,i))
#labels = np.arange(0,15,1)

#element to plot, BASE 0
#style = sys.argv[1] #k (k grid), rc (reaction coordinate), o (other)


filenames = sys.argv[1:]

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')

linestyles = ['-','--']*50

ytop = 1
for i,filename in enumerate(filenames):
    print(filename)
    bins,re,im,dimension,max_e = read_memory_kernel(filename,treat_complex=False)
    print(dimension)
    c=0
    el_idx=[]
    for ii in range(dimension):
        for jj in range(ii,dimension):
            if ii==jj:
                el_idx.append(c)
            c+=1
    print(el_idx)
    color_idx = np.linspace(0, 1, len(el_idx))

    idx = (np.argwhere(bins==1)[0])[0]
    ymax = np.max(re[:,:idx])
    if ymax > ytop:
        ytop = ymax
    output_dir = os.path.dirname(filename)

    for en,d in enumerate(el_idx):
        label = str(en)
        ax.plot(bins,re[d,:],linestyle=linestyles[i],linewidth=0.7,label=label,color=plt.cm.copper(color_idx[len(el_idx)-en-1]))


   

    





#ax.set_yticks(np.arange(0, 2.5, 0.5))


ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.2))


#ax.set_ylim(bottom=0,top=1.)
#ax.set_xlim(0,np.max(bins))

ytop = np.ceil(ytop)

ax.set_ylim(0,ytop)
if ytop < 5:
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
if ytop > 5 and ytop < 15:
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))




ax.set_xlim(left=0,right=1)
#if len(filenames) > 1:
#    ax.legend(loc=1,ncol=1,fancybox=True,framealpha=0)
ax.set_xlabel("Excitation energy / eV")
ax.set_ylabel(r'$\Lambda_{ij}(\epsilon)\ /\ \mathrm{ps}^{-1} $')
fig.set_figheight(2.0)
fig.set_figwidth(3.25)
#plt.gcf().subplots_adjust(left=0.2,bottom=0.2)
#fig.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)
plt.legend(ncol=5,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.1), loc='center')
fig.savefig('memory_elements.pdf',transparent=True,bbox_inches='tight')




