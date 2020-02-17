
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

points=0
with open(sys.argv[1]) as f:
    for line in enumerate(f):
        points+=1

points = points/2
c = 0
coordinate = np.zeros((points))
energy = np.zeros((points))
with open(sys.argv[1]) as f:
    for i,line in enumerate(f):
        if 'aims' in line:
            coordinate[c] = int(line.split()[1])
        elif 'eV' in line:
            energy[c]=float(line.split()[-2])
            c +=1 


fig,ax = plt.subplots(1,1)#sharex=True, facecolor='w')
order = np.argsort(coordinate)
x = coordinate[order]
abs_energy = energy[order]

print(np.amin(abs_energy))
y = np.subtract(abs_energy,np.amin(abs_energy))

ea = np.amax(abs_energy)-np.amin(abs_energy)

ax.plot(x,y,'-',marker='o',mfc='none')
ax.set_xlim(0,points-1)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))

fig.text(0.5,0.5,"Barrier = {0:.2f} eV".format(ea),fontsize=12)
fig.text(0.5, 0.00, r"Reaction coordinate", ha='center',fontsize=15)
fig.text(0.02, 0.5, r'Relative energy / eV', va='center', rotation='vertical',fontsize=15)
fig.savefig('barrier.pdf',transparent=True,bbox_inches='tight')

