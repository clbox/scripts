#!/usr/bin/env python
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



color_idx = np.linspace(0, 1, len(sys.argv[1:]))
fig = plt.figure()
ax = fig.gca(projection='3d')
x_max = 1
bins[bins>1] = np.nan
for i in range(len(time_axis)):
    y = np.ones((len(bins)))*time_axis[i]

    ax.plot(bins,y,time_re[i,0,:], color = plt.cm.cool(color_idx[i]),zorder=len(time_axis)-i)


# # set the axis limits
#ax.set_ylim(0,95)
ax.set_xlim(0,x_max)
# add axis labels
ax.set_xlabel("Excitation energy / eV")
ax.set_ylabel('Time / fs')
ax.set_zlabel(r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $')

fig.savefig('time_3d_projected_kernel.pdf',transparent=True,bbox_inches='tight')
plt.show()

