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

#FOURIER TRANSFOR TO ETA(T)
from ase.units import _hbar, J, s, fs
hbar = _hbar *J #eV s

frequencies=(bins/hbar)/s #ase units inverse time
times = 1/frequencies #ase units time
eta_t = np.zeros((len(sys.argv[1:]),dimension,dimension,len(times)))
for p,path in enumerate(sys.argv[1:]):  
    i_cart = -1
    i_atom = 0
    for i in range(dimension):
        i_cart = i_cart + 1
        if (i_cart>2):
            i_cart = 0
            i_atom = i_atom + 1
                
        j_atom = 0
        j_cart = -1          
        for j in range(dimension):
            j_cart = j_cart + 1
            if (j_cart>2):
                j_cart = 0
                j_atom = j_atom + 1
                    
            if j<i:
                continue     
    
            lambda_omega = re_memory_kernel[p,i,j,:]/(fs*1000) 
            func = np.zeros(len(frequencies))
            for t in range(len(times)):
                for w in range(len(frequencies)):
                    func[w]=(lambda_omega[w])*np.cos(frequencies[w]*times[t])
                        
                    #ase(time)^-2
                eta_t[p,i,j,t]=(1/(2*np.pi))*np.trapz(func[:],frequencies[:])





fig6, ax = plt.subplots(dimension, dimension, sharex='all', sharey='all')

for p,path in enumerate(sys.argv[1:]): 
    for i in range(dimension):
        for j in range(dimension):
            if j>=i:
                ax[i,j].plot(times/fs,eta_t[p,i,j,:],'.',label=str(path),markersize=3,linestyle='-',linewidth=0.1)
                ax[i,j].set_xlim(0,15)

                #index = int(np.where(bins==1)[0])
                #ax[i,j].set_ylim(np.min((re_memory_kernel[:,:,:,:-index])),np.max(re_memory_kernel[:,:,:,:-index]))
    
ax[0,-1].legend()
fig6.set_figheight(12)
fig6.set_figwidth(12)
fig6.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig6.text(0.01, 0.5, r'$\eta(t)\ /\ \AA^{-2} \mathrm{\frac{ev}{amu}} $', va='center', rotation='vertical',fontsize=15)
fig6.savefig('eta_t_tensor.pdf',transparent=True,bbox_inches='tight')
