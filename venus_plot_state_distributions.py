import numpy as np
import sys
import os
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)

filenames = sys.argv[1:]

x3_exp = np.arange(0,7,1)
v3_exp = [0.0,0.2195652173913044,0.42391304347826086,0.35434782608695653,0.0,0,0]

x11_exp = np.arange(2,12)
v11_exp = [0.047,0.099,0.18,0.18,0.18,0.12,0.068,0.052,0.036,0.025]

x15_exp = np.arange(5,16)
v15_exp = [0.115,0.1339,0.194,0.192,0.125,0.082,0.04,0.05,0.019,0.015,0.036]

x16_exp = np.arange(0,17,1)
v16_exp = [0.0,0.0,0.04,0.08,0.13,0.15,0.19,0.11,0.12,0.07,0.04,0.02,0.03,0.02,0.01,0.02,0.02]

tdpt_args = {'marker' : '.','linestyle' : '-','color' : 'purple', 'label' : r'EANN-EFT', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD PES', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-','color' : 'blue', 'label' : r'EANN-LDFA', 'alpha' : 1.0}



plotted_exp = False
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
for i,filename in enumerate(filenames):
    if 'v03' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x3_exp,v3_exp,color='black',label=r'$\nu_i=3$ exp')
        plotted_exp = True

    if 'v11' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x11_exp,v11_exp,color='black',label=r'$\nu_i=11$ exp')
        plotted_exp = True

    if 'v15' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x15_exp,v15_exp,color='black',label=r'$\nu_i=15$ exp')
        plotted_exp = True

    if 'v16' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x16_exp,v16_exp,color='black',label=r'$\nu_i=16$ exp')
        plotted_exp = True

    dis = np.loadtxt(filename)

    if 'tdpt' in os.path.abspath(filename):
        mode_args = tdpt_args

    if 'bomd' in os.path.abspath(filename):
        mode_args = bomd_args

    if 'ldfa' in os.path.abspath(filename):
        mode_args = ldfa_args

    if 'i2' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 2$'
        mode_args['linestyle'] = '--'

    if 'i3' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 2$'
        mode_args['linestyle'] = '--'
    
    if 'i4' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 2$'
        mode_args['linestyle'] = '-.'
    


    a = ax.plot(dis[:,0],dis[:,1],**mode_args)



    #if 'tdpt' in os.path.abspath(filename):
    #    '.-',color='purple',label=r'TDPT')# $\times 2$')
    #marker='^',linestyle='-',color='red',label=r'BOMD')
    # marker='s',linestyle='-',color='blue',label=r'LDFA $\times 4$'))

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.legend()


fig.set_figheight(4)
fig.set_figwidth(5)
fig.text(0.5, 0.00, r"Final vibrational state ($\nu_f$)", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Population', va='center', rotation='vertical',fontsize=15)
fig.savefig('state2state.pdf',transparent=True,bbox_inches='tight')