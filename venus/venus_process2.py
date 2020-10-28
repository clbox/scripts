import numpy as np
from ase.units import _hbar, J, s, fs
from ase import Atoms
import os
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)
import sys
import scipy
from scipy import optimize
import seaborn as sns


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

matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"

filenames = sys.argv[1:]

nstates = len(filenames)

traj_dict = {'mode' : [],
            'angle' : [],
            'lifetime' : [],
            'jf' : [],
            'vf' : [],
            'i_v' : [],
            'i_r' : [],
            'i_t' : [],
            'f_v' : [],
            'f_r' : [],
            'f_t' : [],
                            }

ntrajs=0
mode = ''
ntrapped = 0
for i,filename in enumerate(filenames):
    if 'tdpt' in os.path.abspath(filename):
        mode = 'ODF'
    if 'bomd' in os.path.abspath(filename):
        mode = 'BOMD'
    if 'ldfa' in os.path.abspath(filename):
        mode = 'LDFA'
    with open(filename) as f:
        if 'trapped' in filename:
            vf = -1
            ntrapped+=1
            continue
        else:
            vf = filename.split('/')[-1]
            vf = vf.replace('.dat','')
            vf = int(vf)
        for line in f:
            if 'Trajectory' in line:
                ntrajs += 1
                traj_dict['vf'].append(vf)
                traj_dict['mode'].append(mode)
            elif 'Lifetime' in line:
                numbers = line.replace(',','')
                traj_dict['lifetime'].append(float(numbers.split()[2]))
                traj_dict['angle'].append(float(numbers.split()[7]))
                traj_dict['jf'].append(int(numbers.split()[-1]))
            elif 'Initial' in line:
                numbers = line.replace(',',' ')
                traj_dict['i_v'].append(float(numbers.split()[8]))
                traj_dict['i_r'].append(float(numbers.split()[9]))
                traj_dict['i_t'].append(float(numbers.split()[10]))
            elif 'Final' in line:
                numbers = line.replace(',',' ')
                traj_dict['f_v'].append(float(numbers.split()[8]))
                traj_dict['f_r'].append(float(numbers.split()[9]))
                traj_dict['f_t'].append(float(numbers.split()[10]))


unique = list(set(traj_dict['mode']))
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
for mode in unique:
    idx = (np.where((np.array(traj_dict['mode'])==mode)))[0]
    f_v = np.array(np.array(traj_dict['f_v'])[idx])
    sns.distplot(f_v-0.113,ax=ax,label=mode)
ax.legend()
ax.set_xlim(-0.12,0.8)
ax.set_xlabel(r'$\Delta E_{\mathrm{vib}}$ / eV')
ax.set_yscale('log')
fig.savefig('f_v.pdf',transparent=True,bbox_inches='tight')


fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
for mode in unique:
    idx = (np.where((np.array(traj_dict['mode'])==mode)))[0]
    f_t = np.array(np.array(traj_dict['f_t'])[idx])
    sns.distplot(f_t,ax=ax,label=mode)
ax.legend()
#ax.set_xlim(-0.12,0.8)
ax.set_xlabel(r'$E_{f}$ / eV')
ax.set_yscale('log')
fig.savefig('f_t.pdf',transparent=True,bbox_inches='tight')

fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
cmaps = ['Blues','Reds','Greens']
#for i,mode in enumerate(unique):
for i,mode in enumerate([0,1,2]):
    idx = (np.where((np.array(traj_dict['vf'])==mode)))[0]
    f_r = np.array(np.array(traj_dict['f_r'])[idx])
    f_t = np.array(np.array(traj_dict['f_t'])[idx])
    sns.kdeplot(f_t,f_r,ax=ax,label=mode,cmap=cmaps[i],gridsize=60)
ax.legend()
#ax.set_ylim(0,0.2)
#ax.set_xlim(-0.12,0.8)
ax.set_xlabel(r'$E_{f}$ / eV')
ax.set_xlabel(r'$E_{rot}$ / eV')
#ax.set_yscale('log')
fig.savefig('f_r.pdf',transparent=True,bbox_inches='tight')