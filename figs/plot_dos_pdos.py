import numpy as np
import os
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)

tdos_args = { 'linestyle' : '-','color' : 'black', 'label' : r'Total'}
largs = { 'linestyle' : '-'}
colours = ['dodgerblue','lightgreen','mediumvioletred','orangered']
#Total dos is argument 1
dos_file = sys.argv[1]

#pdos files are rest of arguments if givent
pdos_files = sys.argv[2:]


total_dos = np.loadtxt(dos_file,skiprows=3)

fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)

ldos = np.zeros((len(total_dos),4))

for pdos_file in pdos_files:
    data = np.loadtxt(pdos_file,skiprows=4)

    for column in range(np.shape(data[:,2:])[1]):
        ldos[:,column] += data[:,2+column]




ax.plot(total_dos[:,0],total_dos[:,1],**tdos_args)

for l in range(3):
    ax.plot(total_dos[:,0],ldos[:,l],**largs,color=colours[l],label=r'$l = $'+str(l))

ax.legend(loc=1)
font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)


ax.tick_params(axis='both', which='major', labelsize=12)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

ax.xaxis.set_minor_locator(MultipleLocator(0.5))
ax.xaxis.set_major_locator(MultipleLocator(2.0))

ax.set_xlim(np.min(total_dos[:,0]),np.max(total_dos[:,0]))
ax.set_ylim(bottom=0)
fig.set_figheight(3.0)
fig.set_figwidth(4.25)
#fig.set_constrained_layout_pads(w_pad=0, h_pad=0)

if 'raw' not in dos_file:
    ax.set_xlabel('$E - E_f$ / eV',fontsize=12,fontname=font,color='black')
else:
    ax.set_xlabel('E / eV',fontsize=12,fontname=font,color='black')
ax.set_ylabel('DOS',fontsize=12,fontname=font,color='black')
plt.gcf().subplots_adjust(left=0.2,bottom=0.3)
#fig.text(0.5, 0.00, r"Final vibrational state ($\nu_f$)", ha='center',fontsize=15)
#fig.text(0.01, 0.5, 'Population', va='center', rotation='vertical',fontsize=15)
fig.savefig('dos.pdf',transparent=True)#,bbox_inches='tight')