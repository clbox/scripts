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

energy_files = sys.argv[1:]

linestyles = ['-',':','-.','--']
markers = ['o','^','s','v']
#labels = [r'$E(S = 0)$',r'$E(S = 1)$']
labels = ['0','1','2','3']

height_dirs = np.arange(0,15,1)

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')

lowest_energy=0
with open(energy_files[1]) as f:
     lines = f.readlines()
     y = [float(line.split()[1]) for line in lines]
     y = np.array(y)
     if np.min(y) < lowest_energy:
         lowest_energy = np.min(y)


for i,filename in enumerate(energy_files):
    with open(filename) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]
        x = np.array(x)
        y = np.array(y) - lowest_energy

        if 'eqm' in filename:
            colours = ['blue','dodgerblue']
            labels = [r'$r = 1.3, S=0$',r'$r = 1.3, S=1$']
        else:
            colours = ['red','salmon']
            labels = [r'$r = 1.7, S=0$',r'$r = 1.7, S=1$']
        if 'none' in filename:
            linestyle = '-'
            colour = colours[0]
            label =labels[0]
        else:
            linestyle = '--'
            colour = colours[1]
            label =labels[1]
        #normalise_y = y[-1]
        #y = y-normalise_y
        # y = y - lowest_energy
        ax.plot(x,y,linestyle=linestyle,label=label,color=colour,linewidth=2,marker=markers[i])

#PBE
#ax.set_ylim( -19288900, -19288600.7)
#BEEF
#ax.set_ylim(-19298636,-19298630)
#NORMALISE
ax.set_ylim(-0.1,5.5)

ax.legend()
fig.text(0.5, 0.01, r'COM height above the surface / $\AA$', ha='center',fontsize=15)
fig.text(0.01, 0.5, r"$\mathrm{E_{ads}}$ / eV", va='center',rotation='vertical',fontsize=15)
fig.set_figwidth(6)
fig.set_figheight(5)
fig.savefig('energies.pdf',transparent=True,bbox_inches='tight')







