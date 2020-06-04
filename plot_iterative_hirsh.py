import numpy as np
import sys
import os
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
#Reads in several hirsh.txt and plots


paths = sys.argv[1:]

x = np.linspace(0.1,7.1,15)
x = x[2:]

sum_charges = []
a1_charges = []
a2_charges = []
a3_charges = []

for path in paths:
    charges = np.loadtxt(path)

    sum_charges.append(np.sum(charges))
    a1_charges.append(charges[0])
    a2_charges.append(charges[1])
    a3_charges.append(np.sum(charges[2:]))




fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')


ax.plot(x,a1_charges,color='red',linewidth=0.5,label='Oxygen')
ax.plot(x,a2_charges,color='blue',linewidth=0.5,label='Nitrogen')
ax.plot(x,a3_charges,color='gold',linewidth=0.5,label='Gold')
ax.plot(x,sum_charges,color='black',linewidth=0.5,label='Sum')

ax.legend(loc=1)

ax.set_xlabel(r"Height / $\AA{}$",fontsize=15)
ax.set_ylabel(r'Hirshfeld charge / $e^{-}$',fontsize=15)

fig.set_figheight(4)
fig.set_figwidth(5)
#fig.text(0.5, 0.00, r"Height / $\AA{}$", ha='center',fontsize=15)
#fig.text(0.01, 0.5, r'Hirshfeld charge / $e^{-}$', va='center', rotation='vertical',fontsize=15)
fig.savefig('hirsh.pdf',transparent=True,bbox_inches='tight')
