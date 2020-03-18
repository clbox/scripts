#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys


from mem_energy_loss import read_memory_kernel

mem_file = (glob.glob(sys.argv[1]+'/*friction_memory_kernel.out'))[0]
bins,re,im,dimension,max_e = read_memory_kernel(mem_file)



fig6, ax = plt.subplots(dimension, dimension, sharex='all', sharey='all')
max_value = 0
min_value = 0
for ii,path in enumerate(sys.argv[1:]):
    bins,re,im,dimension,max_e = read_memory_kernel(path)

    if bins[-1] > 1.0:
        index = int(np.where(bins==1)[0])
    else:
        index = -1
    if np.max(re_memory_kernel[:,:,:index]) > max_value:
        max_value = np.max(re_memory_kernel[:,:,:index])
    if np.min((re_memory_kernel[:,:,:index])) < min_value:
        min_value = np.min((re_memory_kernel[:,:,:index]))

    for i in range(dimension):
        for j in range(dimension):
            if j>=i:
                ax[i,j].plot(bins,re_memory_kernel[i,j,:],'.',markersize=3,
                        linestyle='-',linewidth=0.1, label=str(path),
                        color = plt.cm.cool(color_idx[ii]))
                ax[i,j].set_xlim(0,1)
                ax[i,j].set_ylim(min_value,max_value)
                ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.05))

if len(sys.argv[1:]) > 10:
    print('no legend as too many plots')
else:
    ax[0,-1].legend()
fig6.set_figheight(20)
fig6.set_figwidth(20)
fig6.text(0.5, 0.01, "Excitation energy / eV", ha='center',fontsize=15)
fig6.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)
fig6.savefig('memory_kernel_tensor.png',dpi=300,transparent=True,bbox_inches='tight')
fig6.savefig('memory_kernel_tensor.pdf',transparent=True,bbox_inches='tight')