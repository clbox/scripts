import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
import math
import sys
hbar = _hbar * J * s 
ps = fs*1000

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import memory_eta_analysis as mea

line_width = 0.4
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.6
matplotlib.rcParams['lines.linewidth'] = 0.6

masses = mea.get_friction_masses('geometry.in',[16,17])

path = 'friction_gamma2.out'

mem_cutoff = 15 #fs, maximum memcutoff
cutoffs = np.arange(1,mem_cutoff,1) #range of cutoffs to test integral convergence


time_step = 2 #fs not really important, discretised to time_step / 40

fig, ax = plt.subplots(1,1)

integrals = np.zeros(len(cutoffs))

times,eta_t,dimension = mea.get_eta_t(path,mem_cutoff,time_step,masses)

e=0
for i in range(dimension):
    for j in range(i,dimension):
        i_atom,j_atom = i // 3, j // 3
        if i == j == 0:
            ax.plot(times/fs,eta_t[e,:],color='black')
            
            for co,cutoff in enumerate(cutoffs):
                integral = np.trapz(eta_t[e,times<=cutoff*fs],times[times<=cutoff*fs])
                integral /= np.sqrt(masses[i_atom]*masses[j_atom])
                integral*=ps
                integrals[co]=integral
    e+=1

#ax.legend()
ax.set_xlim(0,mem_cutoff)

ax.set_xlabel('Time / fs')
ax.set_ylabel(r'$\Lambda(t)$ / eV $\AA^{-2}$')


ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')


fig.set_figheight(3.5)
fig.set_figwidth(5.)

fig.savefig('eta_t.pdf')