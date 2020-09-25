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
# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
annotate=True
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"



boundaries = True
i = 0
j = 0
labels = ['Adsorption','TS','Product']
colours = ['red','navy','black','orange','orange','orange']
markers = ['o','^','s','^']
filenames = ['adsorption','TS','dissociation']




fig, ax = plt.subplots(1, 2,sharey='all')#, constrained_layout=True)

for m,mode in enumerate(['kgrid','broadening']):
    for c,filename in enumerate(filenames):
        data = np.loadtxt(filename+'/conventional/'+mode+'/projected_tensors.out')

        print(np.shape(data))
        ndim = np.shape(data)[1]

        ntensors = np.shape(data)[0] // ndim

        data = data.reshape(ntensors,ndim,ndim)

        x_data = []
        with open(filename+'/conventional/'+mode+'/projected_tensors.out','r') as f:
            for line in f:
                if '#' in line:
                    line = line.replace('#','')
                    line = line.replace('/',' ')
                    if mode=='kgrid':
                        x_data.append(int(line.split()[0]))
                    else:
                        x_data.append(float(line.split()[0]))

        if c == 0:
            xmin = min(x_data)
            xmax = max(x_data)
        else:
            if min(x_data) < xmin:
                xmin = min(x_data)

            if max(x_data) > xmax:
                xmax = max(x_data)
        

        if boundaries:
            if mode=='kgrid':
                final = data[-1,i,j]
            if mode=='broadening':
                x_data = np.array(x_data)
                final = data[x_data==0.6,i,j]
            ax[m].axhline(y=(final+0.1*final), xmin=0, xmax=100,color=colours[c],linestyle='--',zorder=0)
            ax[m].axhline(y=(final-0.1*final), xmin=0, xmax=100,color=colours[c],linestyle='--',zorder=0)
        ax[m].plot(x_data,data[:,i,j],label=labels[c],color=colours[c],marker=markers[c],mfc='none')



#######################################


#ax.legend(fontsize=12)
annotate_args = {'xy' : (0.03,0.91), 'xycoords' : 'axes fraction'}
ax[0].annotate('(a)',ha="left", **annotate_args)
ax[1].annotate('(b)',ha="left", **annotate_args)

ax[0].yaxis.set_minor_locator(MultipleLocator(0.05))
ax[0].yaxis.set_major_locator(MultipleLocator(0.1))
ax[0].set_ylim(bottom=0.1,top=0.8)
ax[0].set_xlim(8,16)
ax[0].xaxis.set_minor_locator(MultipleLocator(1))

ax[1].set_xlim(0,1)
ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))
ax[1].xaxis.set_major_locator(MultipleLocator(0.2))




fig.set_figheight(2.0)
fig.set_figwidth(6.5)


plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(-0.1, 1.2), loc='center')
#plt.tight_layout()
plt.subplots_adjust(hspace=0.45)



ax[0].set_ylabel(r'$\Lambda_{ij}$ / ps$^{-1}$')#,color='white')
ax[0].set_xlabel('$N_k$')#,color='white')
ax[1].set_xlabel(r'$\sigma$ / eV')#,color='white')
fig.savefig('figs3.pdf',transparent=True,bbox_inches='tight')
fig.savefig('figs3.tiff',transparent=True,bbox_inches='tight',dpi=600)
fig.savefig('figs3.eps',transparent=True,bbox_inches='tight')


