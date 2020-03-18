#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy.interpolate import griddata
import sys
import glob
from mem_energy_loss import read_memory_kernel

mem_file = (glob.glob(sys.argv[1]))[0]
bins,re,im,dimension,max_e = read_memory_kernel(mem_file)

time_step = 1
n_files = len(sys.argv[1:])
print('number of files: ' + str(n_files))
time_axis = np.zeros((n_files))

time_re = np.zeros((len(sys.argv[1:]),np.shape(re)))
for ii,path in enumerate(sys.argv[1:]):
    bins,re,im,dimension,max_e = read_memory_kernel(path)
    time_re[ii,:,:] = re
    time_axis[ii] = float(path.split('/')[0])*time_step

fig, ax = plt.subplots(dimension, dimension, sharex='all', sharey='all')

grid_x, grid_y = np.mgrid[0:1:100j, 0:np.max(time_axis):n_files*1j]


for i in range(dimension):
    for j in range(i,dimension):
        x = bins
        #y = 
        grid_z = griddata(points, values, (grid_x, grid_y), method='nearest')
        


