import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys

for path in sys.argv[1:]:
    head_count =0
    header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction"] #skip lines
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                dimension = int(line.split()[3])
                head_count += 1
            if "Discretization" in line:
                discretization=float(line.split()[-1])
            if any(x in line for x in header):
                continue
            max_e = float(line.split()[0])
    print("Friction max energy = "+str(max_e))
    print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
    elements = (((dimension*dimension)-dimension)/2)+dimension
    print("There are " + str(elements) + " coupling components")
    if elements < head_count:
        n_spin = 2
        print("This system is spin unrestricted")

    bins=np.zeros((int(max_e/discretization)+1))
    print(len(bins))
    re_memory_kernel=np.zeros((len(sys.argv[1:]),dimension,dimension,len(bins)))
    im_memory_kernel=np.zeros_like(re_memory_kernel)
    
for p,path in enumerate(sys.argv[1:]):  
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                i = int(line.split()[3])
                j = int(line.split()[4])
                head_count += 1
                c=0
            if any(x in line for x in header):
                continue
            else:
                re_memory_kernel[p,i-1,j-1,c]=float(line.split()[1])
                im_memory_kernel[p,i-1,j-1,c]=float(line.split()[2])
                bins[c]=float(line.split()[0])
                c +=1

fig6, ax = plt.subplots(dimension, dimension, sharex='all', sharey='all')

for p,path in enumerate(sys.argv[1:]): 
    for i in range(dimension):
        for j in range(dimension):
            if j>=i:
                ax[i,j].plot(bins,re_memory_kernel[p,i,j,:],'.',markersize=3,
                        linestyle='-',linewidth=0.1, label=(str(path)))
                ax[i,j].set_xlim(0,1)
                index = int(np.where(bins==1)[0])
                ax[i,j].set_ylim(np.min((re_memory_kernel[:,:,:,:index])),np.max(re_memory_kernel[:,:,:,:index]))
    
ax[0,-1].legend()
fig6.set_figheight(7)
fig6.set_figwidth(7)
fig6.text(0.5, 0.01, "Excitation energy / eV", ha='center',fontsize=15)
fig6.text(0.01, 0.5, r'$\Lambda(\epsilon)\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)
fig6.savefig('memory_kernel_tensor.pdf',transparent=True,bbox_inches='tight')
