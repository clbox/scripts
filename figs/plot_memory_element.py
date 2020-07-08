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

colours = ['navy','firebrick','dodgerblue','violet']
#labels = np.arange(0,15,1)

#element to plot, BASE 0
element = int(sys.argv[1])

filenames = sys.argv[2:]




fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')

#plot gaussian
x0 = 0
s = 0.6
x = np.linspace(0,6,5000)
c1 = 'maroon'
gauss = (np.exp(-0.5*((x-x0)*(x-x0))/(s*s))/(s*np.sqrt(np.pi)))*(1/np.sqrt(2))
ax.plot(x,gauss,'-',color=c1,linewidth=2)
ax.fill_between(x, gauss,np.zeros_like(gauss),color=c1,alpha=0.4)



for i,filename in enumerate(filenames):
    bins,re,im,dimension,max_e = read_memory_kernel(filename,treat_complex=False)

    c=0
    for ii in range(dimension):
        for jj in range(ii):
            if ii==jj==element:
                break
            c+=1

    ax.plot(bins,re[c,:],linestyle='-',linewidth=1,label=str(labels[i]),color=colours[i])


    output_dir = os.path.dirname(filename)
    if 'project' in filename:
        tensor_file = (glob.glob(output_dir+'/*projected_tensor.out'))[0]
    else:
        tensor_file = (glob.glob(output_dir+'/*friction_tensor.out'))[0]
    tensor = np.loadtxt(tensor_file)
    element_val = tensor[element,element]

    ax.annotate('', xy=(0,element_val), xycoords='data', xytext=(-0.4, element_val), 
            arrowprops=dict(arrowstyle="-|>, head_width=0.3, head_length=0.7",color='navy'),
            #arrowprops=dict(width=0.5),
            color='red')

    c2 = 'black'
    ax.axhline(y=(element_val*4), xmin=0, xmax=100,color=c2,linestyle='--',linewidth=1.5)

    ax.annotate('', xy=(2,4*element_val), xycoords='data', xytext=(2,element_val), 
            arrowprops=dict(arrowstyle="-|>, head_width=0.3, head_length=0.7",color=c2),
            #arrowprops=dict(width=0.5),
            color='red')

    ax.text( x=2.1,y=2*element_val, s=r'$\times 4$', color=c2)
    





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
fig.set_figheight(2.7)
fig.set_figwidth(3.25)
plt.gcf().subplots_adjust(left=0.2,bottom=0.2)
#fig.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)

fig.savefig('memory_element_'+str(element)+'.pdf',transparent=True,bbox_inches='tight')




