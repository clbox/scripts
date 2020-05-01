import numpy as np
import glob
import sys
import os
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


energy_files = ['./none/energies.out','./collinear/energies.out']

linestyles = ['-','--']

labels = ['none','collinear - fixed 1']

dirs = ['none','collinear']

height_dirs = np.arange(0,15,1)

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')
for i,filename in enumerate(energy_files):
    with open(filename) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]
        x = np.array(x)
        x = (x*0.5) + 0.1
        ax.plot(x,y,marker='.',linestyle=linestyles[i],label=labels[i],color='gray',linewidth=2)

ax.set_ylim( -19287077.6, -19287077.1)

ax2 = ax.twinx()
for i,dir in enumerate(dirs):

    heights = []
    tensors = []
    
    spin = str(i)

    for ii,height_dir in enumerate(height_dirs):

        try:
            ft = np.loadtxt(dir+'/be/{:02d}/friction_tensor.out'.format(height_dir))
        except:
            continue
        
        heights.append((height_dir*0.5)+0.1)
        tensors.append(ft)


    heights = np.array(heights)
    tensors = np.array(tensors)
    
    # ax2.plot(heights,tensors[:,0,0],linestyle=linestyles[i],label=r'$\Lambda_{O_x O_x}(S=$'+spin+r'$)$',color='blue')
    # ax2.plot(heights,tensors[:,1,1],linestyle=linestyles[i],label=r'$\Lambda_{O_y O_y}(S=$'+spin+r'$)$',color='red')
    # ax2.plot(heights,tensors[:,2,2],linestyle=linestyles[i],label=r'$\Lambda_{O_z O_z}(S=$'+spin+r'$)$',color='green')
    ax2.plot(heights,tensors[:,0,0],linestyle=linestyles[i],label=r'$\Lambda_{N_x N_x}(S=$'+spin+r'$)$',color='deepskyblue')
    ax2.plot(heights,tensors[:,1,1],linestyle=linestyles[i],label=r'$\Lambda_{N_y N_y}(S=$'+spin+r'$)$',color='orangered')
    ax2.plot(heights,tensors[:,5,5],linestyle=linestyles[i],label=r'$\Lambda_{N_z N_z}(S=$'+spin+r'$)$',color='palegreen')



ax2.set_ylim(0,1.6)
ax.set_xlim(0,7.1)

ax.legend(loc=2)
ax2.legend(loc=0)
fig.text(0.5, 0.01, r'Height above the surface / $\AA$', ha='center',fontsize=15)
fig.text(0.01, 0.5, "Energy / eV", va='center',rotation='vertical',fontsize=15)
fig.text(1.0, 0.5, r"$\Lambda_{ij}$ / ps$^{-1}$", va='center',rotation='vertical',fontsize=15)
fig.set_figwidth(6)
fig.set_figheight(6)
fig.savefig('energies.pdf',transparent=True,bbox_inches='tight')







