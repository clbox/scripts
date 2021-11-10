import gamma.friction_coupling as fc
import sys
import numpy as np
import glob



aims1 = 'aims.out'
gamma_files1 = glob.glob("*gamma*k*out")
gamma_files1.sort()

a = fc.friction_gamma_parser(aims_file=aims1,gamma_files=gamma_files1)

sigmas = np.linspace(0.0001,1,50)

modes = ['default','double_delta','prb_print','no_norm']

for mode in modes:
    tensors = []
    for sigma in sigmas:
        b = fc.friction_tensor(a,0.0,sigma,nspin=1)
        tensor = b.calc_tensor(mode=mode)
        tensors.append(tensor)
    
        dimension = np.shape(tensor)[0]
    
    
    tensors = np.array(tensors)
    for i in range(dimension):
        for j in range(dimension):
            if i==j:
                np.savetxt('gamma_diagonal_'+str(i)+'_'+mode+'.txt',tensors[:,i,j])
