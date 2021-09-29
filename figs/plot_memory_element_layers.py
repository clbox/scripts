#Ideally takes element number as first argument then a list of filenames of *_memory_kernel.out 

#Additionally can read in tensors and mark where the tensor evaluated to.

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
import glob
from mem_energy_loss import read_memory_kernel

line_width = 0.4
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.6
matplotlib.rcParams['lines.linewidth'] = 0.6

#labels = list(range(0,15))

colours = ['lightgrey','grey','black']

labels =[]
# for i in [2,3,4]:
#     labels.append(r'p({}$\times${})'.format(i,i))
#labels = np.arange(0,15,1)

#element to plot, BASE 0
#style = sys.argv[1] #k (k grid), rc (reaction coordinate), o (other)

#element = int(sys.argv[2])

filenames = sys.argv[1:]

filenames.sort()
print(filenames)


n_layers = []
for filename in filenames:
    n_layers.append(int(os.path.dirname(filename)))

max_layers = n_layers[-1]

fig, ax = plt.subplots(6, 1, sharex='all')#, sharey='all')
#color_idx = np.linspace(0, 1, len(n_layers))
linestyles = ['-']*100


ytop = 1
for i,filename in enumerate(filenames):
    print(filename)
    bins,re,im,dimension,max_e = read_memory_kernel(filename,treat_complex=False)
    

    for element in range(6): 


        label = os.path.dirname(filename)
        c=0
        for ii in range(dimension):
            for jj in range(ii,dimension):
                if ii==element and jj==element:
                    print('c' + str(c))
                    d = c
                    break
                c+=1
        print(d)

        idx = (np.argwhere(bins==1)[0])[0]
        ymax = np.max(re[d,:idx])
        if ymax > ytop:
            ytop = ymax
        output_dir = os.path.dirname(filename)

        ax[element].plot(bins,re[d,:],linestyle=linestyles[i],linewidth=0.7,label=label,color=plt.cm.Blues(n_layers[i]/max_layers))


   

    


for element in range(6):
    ax[element].xaxis.set_minor_locator(MultipleLocator(0.1))
    ax[element].xaxis.set_major_locator(MultipleLocator(0.2))

    ytop = np.ceil(ytop)

    ax[element].set_ylim(0,ytop)
    if ytop < 5:
        ax[element].yaxis.set_minor_locator(MultipleLocator(0.1))
        ax[element].yaxis.set_major_locator(MultipleLocator(0.5))
    if ytop > 5 and ytop < 15:
        ax[element].yaxis.set_minor_locator(MultipleLocator(0.5))
        ax[element].yaxis.set_major_locator(MultipleLocator(1.0))



    ax[element].xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
    ax[element].xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
    ax[element].yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
    ax[element].yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

    ax[element].set_xlim(left=0,right=1)
    #if len(filenames) > 1:
    #    ax.legend(loc=1,ncol=1,fancybox=True,framealpha=0)
    ax[element].set_xlabel("Excitation energy / eV")
    ax[element].set_ylabel(r'$\Lambda_{ij}(\epsilon)\ /\ \mathrm{ps}^{-1} $')
    #ax[element].legend(fancybox=True,framealpha=0)


fig.set_figheight(25)
fig.set_figwidth(5)
#plt.gcf().subplots_adjust(left=0.2,bottom=0.2)
#fig.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)
plt.legend(ncol=5,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.1), loc='center')

fig.savefig('memory_element_layers.pdf',transparent=True,bbox_inches='tight')




