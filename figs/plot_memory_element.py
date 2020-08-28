#Ideally takes element number as first argument then a list of filenames of *_memory_kernel.out 

#Additionally can read in tensors and mark where the tensor evaluated to.

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
import glob
from mem_energy_loss import read_memory_kernel


#labels = ['Reactant','Adsorption','Transition state','Dissociated']
#labels = ['Adsorption','Transition state','Dissociated']

#labels = list(range(0,15))

colours = ['navy','firebrick','dodgerblue','violet']
#labels = np.arange(0,15,1)

#element to plot, BASE 0
element = int(sys.argv[1])

filenames = sys.argv[2:]




fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')




for i,filename in enumerate(filenames):
    bins,re,im,dimension,max_e = read_memory_kernel(filename,treat_complex=False)

    c=0
    for ii in range(dimension):
        for jj in range(ii):
            if ii==jj==element:
                break
            c+=1
    output_dir = os.path.dirname(filename)
    ax.plot(bins,re[c,:],linestyle='-',linewidth=1,label=str(output_dir),color=colours[i])


   

    





#ax.set_yticks(np.arange(0, 2.5, 0.5))

font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)


ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.3))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))

ax.set_ylim(bottom=0,top=3.0)
ax.set_xlim(0,np.max(bins))

ax.set_xlim(left=0)
if len(filenames) > 1:
    ax.legend(loc=1,ncol=1,fontsize=12,fancybox=True,framealpha=0)
ax.set_xlabel("Excitation energy / eV",fontsize=12,fontname=font)
ax.set_ylabel(r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $',fontsize=12,fontname=font)
fig.set_figheight(2.0)
fig.set_figwidth(3.25)
plt.gcf().subplots_adjust(left=0.2,bottom=0.2)
#fig.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)

fig.savefig('memory_element_'+str(element)+'.pdf',transparent=True,bbox_inches='tight')




