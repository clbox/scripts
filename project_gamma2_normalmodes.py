import numpy as np
import sys
from ase.io import read
import glob
from mem_energy_loss import read_memory_kernel

Header = """No of components        {}    
Discretization length in eV   {}
Number of Bins      {}  
Excitation energy in eV   Re(Coupling element) 
============================================================================\n"""


mem_file = sys.argv[1]
bins,re,im,ndim,max_e = read_memory_kernel(mem_file)

normalmodes =  open(sys.argv[2]).readlines()


modes = np.zeros([ndim,ndim])
for i, line in enumerate(normalmodes):
  if 'Zero-point energy' in line:
      j = i+1
      for a in range(ndim):
          for b in range(ndim):
              modes[a,b] = \
                      float(normalmodes[j].split()[b])
          j += 1

for i in range(ndim):
    modes[i,:]/=np.linalg.norm(modes[i,:])

print(re)

normal_re = np.zeros_like(re)
for b, bin in enumerate(bins):

    #transform friction matrix into normalmode space
    normal_re[b,:,:] = np.dot(modes,np.dot(re[b,:,:],modes.transpose()))



myfile = open('normalmode_friction_gamma2.out', 'w')
myfile.write(Header.format(ndim,bins[1]-bins[0],len(bins)))
for i in range(ndim):
    for j in range(i,ndim):
        myfile.write('Friction component for {} {} \n'.format(i+1,j+1))
        for e, bin in enumerate(bins):
            myfile.write('{:0.4f}  {:0.6f}\n'.format(bin,normal_re[e,i,j]))
myfile.close()



# A = np.dot(modes,np.dot(friction_tensor,modes.transpose()))


