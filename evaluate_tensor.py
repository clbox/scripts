import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
from ase.io import read
import math
import memory_eta_analysis as mea
import sys
import os
hbar = _hbar * J * s 
ps = fs*1000

#read in xxxx_memory_kernel.out(s) (argument 2) and evaluate the tensor at chosen excitation energy (argument 1)

excit_e = float(sys.argv[1]) #in eV
output_name = 'evaluated_tensor_'+str(excit_e)+'.out'
sigma = 0.6 #eV


def gaussian(x,x0,s):
    try:
        a = 1./np.sqrt(2*np.pi*s*s)
        b = np.exp(-0.5*(((x-x0)/s)**2))
        return a*b
    except OverflowError:
        return 0.0


filenames = sys.argv[2:]
for filename in filenames:

    output_dir = os.path.dirname(filename)

    bins,re_memory_kernel,im_memory_kernel,dimension,max_e = mea.read_memory_kernel(filename,treat_complex=False)

    elements = np.shape(re_memory_kernel)[0]
    tensor = np.zeros((dimension,dimension))

    e=0
    for i in range(dimension):
        for j in range(i,dimension):
            delta = gaussian(bins,excit_e,sigma)
            tensor[i,j] = np.sum(re_memory_kernel[e,:]*delta)/np.sum(delta)
            tensor[j,i] = tensor[i,j]
            e+=1

    np.savetxt(output_dir+output_name,tensor,delimiter='  ',fmt='%.4e')


