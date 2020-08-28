import numpy as np
import sys
from ase.io import read
import glob
from mem_energy_loss import read_memory_kernel
import venus.venus_traj_project as vjp

# mem_file = (glob.glob(sys.argv[1]+'/*friction_memory_kernel.out'))[0]
# bins,re,im,dimension,max_e = read_memory_kernel(mem_file)

#projected_kernel = np.zeros((len(sys.argv[1:]),len(bins),dimension,dimension))

Header = """No of components        {}    
Discretization length in eV   {}
Number of Bins      {}  
Excitation energy in eV   Re(Coupling element) 
============================================================================\n"""

for output_dir in sys.argv[1:]:
    mem_file = (glob.glob(output_dir+'/*friction_gamma2.out'))[0]
    print('-------------'+mem_file+'----------------')
    bins,re,im,dimension,max_e = read_memory_kernel(mem_file)

    geo_file =  (glob.glob(output_dir+'/*geometry.in'))[0]
    print('-------------'+geo_file+'----------------')
    atoms = read(geo_file)

    output = (glob.glob(output_dir+'/*aims.out'))[0]
    print('-------------'+output+'----------------')

    friction_atoms = []
    with open(output, "r") as lines:
        for i,line in enumerate(lines):

            if 'Found friction request for atom' in line:
                friction_atoms.append(int(line.split()[-1])-1)
            elif 'Input structure read successfully' in line:
                break
            else:
                continue
    

    modes = vjp.get_modes(atoms,friction_atoms,mode=2)

    re_tensor = np.zeros((len(bins),dimension,dimension))
    e=0
    for i in range(dimension):
        for j in range(i,dimension):
            re_tensor[:,i,j] = re[e,:]
            re_tensor[:,j,i] = re[e,:]
            e+=1

    projected_kernel = np.zeros((len(bins),dimension,dimension))
    for e,bin in enumerate(bins):
        projected_kernel[e,:,:] = np.dot(modes,np.dot(re_tensor[e,:,:],modes.transpose()))

    myfile = open(output_dir+'/projected_memory_kernel.out', 'w')
    myfile.write(Header.format(dimension,bins[1]-bins[0],len(bins)))
    for i in range(dimension):
        for j in range(i,dimension):
            myfile.write('Friction component for {} {} \n'.format(i+1,j+1))
            for e, bin in enumerate(bins):
                myfile.write('{:0.4f}  {:0.6f}\n'.format(bin,projected_kernel[e,i,j]))
    myfile.close()
    