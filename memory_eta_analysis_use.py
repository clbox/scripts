import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
import math
hbar = _hbar * J * s 
ps = fs*1000

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import memory_eta_analysis as mea

masses = mea.get_friction_masses('geometry.in',[16,17])

path = 'friction_memory_kernel.out'

mem_cutoff = 20 #fs

time_step = 2 #fs not really important, discretised to time_step / 40

times,eta_t,dimension = mea.get_eta_t(path,mem_cutoff,time_step,masses)

fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
e=0
for i in range(dimension):
    for j in range(i,dimension):
        ax[i,j].plot(times/fs,eta_t[e,:])
        e+=1
fig.set_figheight(20)
fig.set_figwidth(20)
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'$\eta(t)$ / eV $\AA^{-2}$', va='center', rotation='vertical',fontsize=15)
fig.savefig('eta_t.pdf')