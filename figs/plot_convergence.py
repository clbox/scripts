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

boundaries = True
i = int(sys.argv[1])
j = int(sys.argv[2])
labels = ['Reactant','Adsorption','Transition state','Dissociation']
colours = ['red','navy','black','orange','orange','orange']
markers = ['.','s','o','^']
filenames = sys.argv[3:]
kgrid = False
broad = False
if 'k_grid' in  os.path.abspath(filenames[1]) or 'kgrid' in  os.path.abspath(filenames[1]):
    kgrid = True
    print('kgrid')
if 'broadening' in  os.path.abspath(filenames[1]):
    broad = True
    print('broac')


fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)


for c,filename in enumerate(filenames):
    data = np.loadtxt(filename)

    print(np.shape(data))
    ndim = np.shape(data)[1]

    ntensors = np.shape(data)[0] // ndim

    data = data.reshape(ntensors,ndim,ndim)

    x_data = []
    with open(filename,'r') as f:
        for line in f:
            if '#' in line:
                line = line.replace('#','')
                line = line.replace('/',' ')
                if kgrid:
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
    
    if c ==0:
        continue
    if boundaries:
        if kgrid:
            final = data[-1,i,j]
        if broad:
            x_data = np.array(x_data)
            final = data[x_data==0.6,i,j]
        ax.axhline(y=(final+0.1*final), xmin=0, xmax=100,color=colours[c],linestyle='--')
        ax.axhline(y=(final-0.1*final), xmin=0, xmax=100,color=colours[c],linestyle='--')
    ax.plot(x_data,data[:,i,j],label=labels[c],color=colours[c],marker=markers[c])



#######################################
font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)

ax.tick_params(axis='both', which='major', labelsize=12)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=12)

if kgrid:
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(bottom=0.3,top=0.6)
else:
     ax.set_xlim(0,1)
     ax.xaxis.set_minor_locator(MultipleLocator(0.05))
     ax.xaxis.set_major_locator(MultipleLocator(0.2))
     ax.set_ylim(bottom=0.3,top=0.6)

# ax.xaxis.set_minor_locator(MultipleLocator(0.05))
# ax.xaxis.set_major_locator(MultipleLocator(0.1))



fig.set_figheight(2.0)
fig.set_figwidth(3.25)

plt.gcf().subplots_adjust(left=0.3,bottom=0.3)



ax.set_ylabel(r'$\Lambda_{ij}$ / ps$^{-1}$',fontsize=12,fontname=font)#,color='white')
if kgrid:
    ax.set_xlabel('$N_k$',fontsize=12,fontname=font)#,color='white')
    fig.savefig('kgrid.pdf',transparent=True)#,bbox_inches='tight')
elif broad:
    ax.set_xlabel(r'$\sigma$ / eV',fontsize=12,fontname=font)#,color='white')
    fig.savefig('broad.pdf',transparent=True)#,bbox_inches='tight')
else:
    fig.savefig('convergence.pdf',transparent=True)#,bbox_inches='tight')


fig.legend(ncol=6,fontsize=12,fancybox=True,framealpha=0)
fig.set_figwidth(7)
ax.remove()
fig.savefig('legend.pdf',transparent=True)#,bbox_inches='tight')