import numpy as np
import sys
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)



filenames = sys.argv[1:]

e_diff_list = np.load(filenames[0],allow_pickle=True)
nstates = np.shape(e_diff_list)[0]

colours = ['purple','red']

#fig2, ax2 = plt.subplots(nstates,1, sharex='all',sharey='all')
fig2, ax2 = plt.subplots(3,1, sharex='all',sharey='all')
for f,filename in enumerate(filenames):
    e_diff_list = np.load(filename,allow_pickle=True)
    c = colours[f]

    #for i in range(nstates):
    for i in range(3):
        ax2[i].set_xlim(0.5,3.5)
        ax2[i].boxplot(e_diff_list[i],showfliers=False,
        boxprops=dict(color=c),whiskerprops=dict(color=c),capprops=dict(color=c),medianprops=dict(color=c,linestyle='--'))
        ax2[i].text(s='Nf = {}'.format(i+1),x=0.5,y=0.8)
        #ax2[i].text(s='ntrajs = '+str(ntraj_list[i]),x=0.5,y=0.5)


plt.xticks([1,2,3],['Vibrational','Rotational','Translational'])
#fig2.set_figheight(10*nstates/10)
fig2.set_figheight(5)
fig2.set_figwidth(4)
fig2.text(0.5, 0.03, "Component", ha='center',fontsize=15)
fig2.text(0.01, 0.5, r'E$_i$ - E$_f$ / eV', va='center', rotation='vertical',fontsize=15)
fig2.savefig('energy_distribution_compare.pdf',transparent=True,bbox_inches='tight')