import numpy as np
import sys
from ase.io import read
import glob


def build_tensor(raw_tensor):
    a = len(raw_tensor[0].split())
    lines_per_line = 0 
    for i,line in enumerate(raw_tensor):
        if i == 0:
            continue
        if len(line.split()) != a:
            lines_per_line += 1
        if len(line.split()) == a:
            break

    ndim = (a - 1)*(lines_per_line+1)
    friction_tensor = np.zeros((ndim,ndim))
    c=0
    i=-1
    for ii, line in enumerate(raw_tensor):
        
        if c == lines_per_line+1 or c==0:
            i+=1
            c=0

        for j in range(int(ndim/(lines_per_line+1))):
            if c ==0:
                friction_tensor[i,j]=float(line.split()[j+1])
            else:
                friction_tensor[i,j+((int(ndim/(lines_per_line+1)))*c)]=float(line.split()[j])
        c+=1

    return friction_tensor,ndim

def get_modes(atoms,friction_atoms,mode=2):

    if mode == 1:
        modes = calc_modes(atoms,friction_atoms)

    elif mode == 2:
        modes = calc_modes2(atoms,friction_atoms)

    return modes
    
def calc_modes(atoms,friction_atoms):
        """Calculates required transformation matrix to convert diatomic
        friction tensor to internal coordinates as defined in Maurer et al PRL 2017
        """
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
        #return modes.transpose() #test

def calc_modes2(atoms,friction_atoms):
        """Calculates required transformation matrix to convert diatomic
        friction tensor to internal coordinates as defined by B. Jiang
        """

        ndim = len(friction_atoms)*3
        modes = np.zeros([ndim,ndim])
        f1 = friction_atoms[0]
        f2 = friction_atoms[1]

        pos1 = atoms.positions[f1]
        pos2 = atoms.positions[f2]

        x1 = pos1[0]
        y1 = pos1[1]
        z1 = pos1[2]
        x2 = pos2[0]
        y2 = pos2[1]
        z2 = pos2[2]

        m1 = atoms.get_masses()[f1]
        m2 = atoms.get_masses()[f2]

        mt = m1 + m2

        mr1 = m1/mt

        mr2 = m2/mt

        r = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        r1 = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        #mode 1 - r
        modes[0,0] = ((x1-x2)*mr2)/r
        modes[1,0] = ((y1-y2)*mr2)/r
        modes[2,0] = ((z1-z2)*mr2)/r

        modes[3,0] = ((x2-x1)*mr1)/r
        modes[4,0] = ((y2-y1)*mr1)/r
        modes[5,0] = ((z2-z1)*mr1)/r

        #mode 2  - theta
        modes[0,1] = ((x1-x2)*(z1-z2)*mr2)/r1
        modes[1,1] = ((y1-y2)*(z1-z2)*mr2)/r1
        modes[2,1] = -r1*mr2

        modes[3,1] = ((x2-x1)*(z1-z2)*mr1)/r1
        modes[4,1] = ((y2-y1)*(z1-z2)*mr1)/r1
        modes[5,1] = r1*mr1

        #mode 3 - phi
        modes[0,2] = -(y1-y2)*mr2
        modes[1,2] = (x1-x2)*mr2
        modes[2,2] = 0

        modes[3,2] = -(y2-y1)*mr1
        modes[4,2] = (x2-x1)*mr1
        modes[5,2] = 0

        #mode 4 is the x translation
        modes[:,3] = [1.,0.,0.,1.,0.,0.]
        #mode 5 is the y translation
        modes[:,4] = [0.,1.,0.,0.,1.,0.]
        #mode 6 is the z translation
        modes[:,5] = [0.,0.,1.,0.,0.,1.]

        return modes.transpose()
        #return modes

for output_dir in sys.argv[1:]:
    output = (glob.glob(output_dir+'/*aims.out'))[0]
    print('-------------'+output+'----------------')
    geo_file =  (glob.glob(output_dir+'/*geometry.in'))[0]
    print('-------------'+geo_file+'----------------')


    friction_atoms = []
    raw_tensor = []
    ptensor = False
    pvecs = False
    raw_vecs = []
    with open(output, "r") as lines:
        for i,line in enumerate(lines):

            if 'Found friction request for atom' in line:
                friction_atoms.append(int(line.split()[-1])-1)

            elif '********Printing Friction Tensor in 1/ps*********' in line:
                ptensor = True
                continue
            elif '**********END Printing Friction Tensor************' in line:
                ptensor = False

            elif ptensor:
                raw_tensor.append(line)

            elif '****Diagonalized Lifetime Eigvecs****' in line:
                pvecs = True

            elif '|------------------------------------' in line:
                pvecs = False

            elif 'Final output of selected total energy values:' in line:
                raw_vecs.pop()
                raw_vecs.pop()
                break
            
            elif pvecs:
                raw_vecs.append(line)

    print('Friction atoms:' + str(friction_atoms))

    friction_tensor,ndim = build_tensor(raw_tensor)
    friction_vecs, ndim = build_tensor(raw_vecs)
    friction_tensor *= 1.66054e-27
    atoms = read(geo_file)
    modes = get_modes(atoms,friction_atoms,mode=1)

    E,V  = np.linalg.eig(friction_tensor[:,:])
    E = sorted(E)
    diag = np.diag(E)
    V = friction_vecs.transpose()
    string = ''
    for e in E:
        string += str(abs(e)) + ' '
    #string += '\n'
    #string += '    d          phi          theta            X             Y             Z\n'
    for i,e in enumerate(E):
        mode = modes[i,:]/np.linalg.norm(modes[i,:])
        string +=  str(((np.dot(mode,np.dot(friction_tensor[:,:],mode))))) + ' '
    print(string)#+'\n'




    #transform friction matrix into normalmode space
    A = np.dot(modes,np.dot(friction_tensor,modes.transpose()))

    print('Diagonal lifetimes')
    for i in range(len(E)):
        print (1./A[i,i])

    print('A_trace')
    print(A.trace())
    print('friction_tensor_trace')
    print(friction_tensor.trace())

    print('Projected_tensor')
    for i in range(ndim):
        string = ''
        for j in range(ndim):
            #string += ' {0:14.8f} '.format(A[i,j])
            string += ' {0:14.8e} '.format(A[i,j])
        print(string)

    np.savetxt(output_dir+'/projected_tensor.out',A,delimiter='  ',fmt='%.4e')