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
labels = ['Adsorption','Transition state','Dissociated']

#labels = list(range(0,15))

colours = ['firebrick','dodgerblue','violet']
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

    ax.plot(bins,re[c,:],linestyle='-',linewidth=0.7,label=str(labels[i]),color=colours[i])


    output_dir = os.path.dirname(filename)
    if 'project' in filename:
        tensor_file = (glob.glob(output_dir+'/*projected_tensor.out'))[0]
    else:
        tensor_file = (glob.glob(output_dir+'/*friction_tensor.out'))[0]
    tensor = np.loadtxt(tensor_file)
    element_val = tensor[element,element]

    ax.annotate('', xy=(0,element_val), xycoords='data', xytext=(-0.2, element_val), 
            arrowprops=dict(arrowstyle="->",color=colours[i]))


#plot gaussian
x0 = 0
s = 0.6
x = bins
gauss = (np.exp(-0.5*((x-x0)*(x-x0))/(s*s))/(s*np.sqrt(np.pi)))*(1/np.sqrt(2))
ax.plot(bins,gauss,'--',color='black')
ax.fill_between(bins, gauss,np.zeros_like(gauss),color='lightgrey',alpha=0.7)


#ax.set_yticks(np.arange(0, 2.5, 0.5))

font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)


ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))

ax.set_ylim(bottom=0,top=5.5)
#ax.set_xlim(0,np.max(bins))

ax.set_xlim(0,2)
ax.legend(loc=1,ncol=1,fontsize=12,fancybox=True,framealpha=0)
ax.set_xlabel("Excitation energy / eV",fontsize=12,fontname=font)
ax.set_ylabel(r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $',fontsize=12,fontname=font)
fig.set_figheight(2.7)
fig.set_figwidth(3.25)
plt.gcf().subplots_adjust(left=0.2,bottom=0.2)
#fig.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)

fig.savefig('memory_element_'+str(element)+'.pdf',transparent=True,bbox_inches='tight')




