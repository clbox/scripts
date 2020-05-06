from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy.interpolate import griddata
import sys
import glob
from mem_energy_loss import read_memory_kernel

mem_file = (glob.glob(sys.argv[1]))[0]
bins,re,im,dimension,max_e = read_memory_kernel(mem_file,treat_complex=False)

time_step = 1
n_files = len(sys.argv[1:])
print('number of files: ' + str(n_files))
time_axis = np.zeros((n_files))

time_re = np.zeros((len(sys.argv[1:]),np.shape(re)[0],np.shape(re)[1]))
for ii,path in enumerate(sys.argv[1:]):
    bins,re,im,dimension,max_e = read_memory_kernel(path,treat_complex=False)
    time_re[ii,:,:] = re
    time_axis[ii] = float(path.split('/')[0])*time_step


xv, yv = np.meshgrid(bins,time_axis)


fig, ax = plt.subplots(dimension, dimension, sharex='all', sharey='all')
e=0
lvls = np.linspace(np.min(time_re),np.max(time_re),20)
for i in range(dimension):
    for j in range(i,dimension):
        im = ax[i,j].contourf(xv,yv,time_re[:,e,:],levels=lvls)
        e+=1

cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)

cbar.set_ticks(np.arange(-0.5, np.max(lvls), 0.5))

 
fig.text(0.5, 0.01, "Excitation energy / eV", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Time / fs', va='center', rotation='vertical',fontsize=15)
fig.text(0.9, 0.5, r'Relaxation rate / $\mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)
ax[0,0].set_xlim(0,1)
fig.savefig('memory_kernel_map.png',dpi=300,bbox_inches='tight')
