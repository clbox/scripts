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
            
            elif pvecs:
                raw_vecs.append(line)



    friction_tensor,ndim = build_tensor(raw_tensor)
    friction_vecs, ndim = build_tensor(raw_vecs)

    modes = np.zeros([ndim,ndim])

    ###now we need to generate modes depending on the H2 geometry
    atoms = read(geo_file)
    #the H atoms are 0 and 1
    pos = atoms.positions

    fa1 = friction_atoms[0]
    fa2 = friction_atoms[1]

    internals = np.zeros(6)
    #com = (pos[0,:] + pos[1,:])/2.
    #vec = pos[0,:] - pos[1,:]
    com = (pos[fa1,:] + pos[fa2,:])/2.
    vec = pos[fa1,:] - pos[fa2,:]


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
            string += ' {0:14.8f} '.format(A[i,j])
        print(string)

    np.savetxt(output_dir+'/projected_tensor.out',A,delimiter='  ')