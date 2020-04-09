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

masses = mea.get_friction_masses('geometry.in',[16,17])

paths = sys.argv[1:]

mem_cutoff = 5 #fs

time_step = 2 #fs not really important, discretised to time_step / 40

fig, ax = plt.subplots(2,1)

integrals = np.zeros((len(paths)))
vals = np.zeros((len(paths)),dtype=int)
for p,path in enumerate(paths):
    times,eta_t,dimension = mea.get_eta_t(path,mem_cutoff,time_step,masses)
    val = int((path.split('/'))[0])
    vals[p] = val
    e=0
    for i in range(dimension):
        for j in range(i,dimension):
            if i == j == 0:
                ax[0].plot(times/fs,eta_t[e,:],label=val)
                integral = np.trapz(eta_t[e,:],times/fs)
                integrals[p]=integral
        e+=1
ax[0].legend()
ax[0].set_xlim(0,mem_cutoff)

fig.text(0.5, 0.5, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.75, r'$\eta(t)$ / eV $\AA^{-2}$', va='center', rotation='vertical',fontsize=15)

indices = np.argsort(vals)
vals = vals[indices]
integrals = integrals[indices]

for p,path in enumerate(paths):
    ax[1].plot(vals,integrals,'-o')


ax[1].set_xlim(np.min(vals),np.max(vals))
fig.set_figheight(10)
fig.set_figwidth(7)
fig.text(0.5, 0.01, r"N${_k}$", ha='center',fontsize=15)
fig.text(0.01, 0.25, r'$\int \eta(t) dt$', va='center', rotation='vertical',fontsize=15)
fig.savefig('eta_t.pdf')