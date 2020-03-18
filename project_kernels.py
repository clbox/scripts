import numpy as np
import sys
from ase.io import read
import glob
from mem_energy_loss import read_memory_kernel

def calc_modes(atoms,friction_atoms):

    ndim = len(friction_atoms)*3
    modes = np.zeros([ndim,ndim])
    internals = np.zeros(6)
    f1 = friction_atoms[0]
    f2 = friction_atoms[1]

    pos1 = atoms.positions[f1]
    pos2 = atoms.positions[f2]

    mass1 = atoms.get_masses()[f1]
    mass2 = atoms.get_masses()[f2]
    
    com = (mass1*pos1[:] + mass2*pos2[:])/(mass1+mass2)



    vec = pos1[:] - pos2[:]
    norm = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]
    internals[0] = np.sqrt(norm)
    #phi
    internals[1] = np.arccos((vec[0]/np.sqrt(vec[0]*vec[0]+vec[1]*vec[1])))
    #theta
    vx_tilde = np.cos(internals[1])*vec[0]+np.sin(internals[1])*vec[1]
    internals[2] = np.arccos((vec[2]/np.sqrt(vec[2]*vec[2]+vx_tilde*vx_tilde)))
    internals[3:] = com[:]
    #mode 1 is the internal stretch
    modes[0,:3] = vec
    modes[0,3:] = -vec
    #mode 2 is 1st rotation angle
    # theta angle between molecular axis and z axis
    modes[1,:] = [-vec[1],vec[0],0.,vec[1],-vec[0],0.]
    #mode 3 is the 2nd rotation angle
    #it is defined as the crossing line between 2 planes
    #plane 1 is defined by the z-axis and the H2 axis
    #plane 2 is defined with the H2 axis as the normal vector
    vec2 = [0.,0.,0.]
    vec2[0] = np.cos(internals[1])*vec[2]
    vec2[1] = np.sin(internals[1])*vec[2]
    vec2[2] = -vx_tilde
    modes[2,:3] = np.array(vec2)
    modes[2,3:] = -np.array(vec2)

    #mode 4 is the x translation
    modes[3,:] = [1.,0.,0.,1.,0.,0.]
    #mode 5 is the y translation
    modes[4,:] = [0.,1.,0.,0.,1.,0.]
    #mode 6 is the z translation
    modes[5,:] = [0.,0.,1.,0.,0.,1.]

    for i in range(ndim):
        modes[i,:]/=np.linalg.norm(modes[i,:])
    return modes


# mem_file = (glob.glob(sys.argv[1]+'/*friction_memory_kernel.out'))[0]
# bins,re,im,dimension,max_e = read_memory_kernel(mem_file)

#projected_kernel = np.zeros((len(sys.argv[1:]),len(bins),dimension,dimension))

Header = """No of components        {}    
Discretization length in eV   {}
Number of Bins      {}  
Excitation energy in eV   Re(Coupling element) 
============================================================================\n"""

for output_dir in sys.argv[1:]:
    mem_file = (glob.glob(output_dir+'/*friction_memory_kernel.out'))[0]
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
    

    modes = calc_modes(atoms,friction_atoms)

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
    