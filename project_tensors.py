import numpy as np
import sys
from ase.io import read
import glob
import venus.venus_traj_project as vjp


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
    try:
        output = (glob.glob(output_dir+'/*aims.out'))[0]
        print('-------------'+output+'----------------')
        geo_file =  (glob.glob(output_dir+'/*geometry.in'))[0]
        print('-------------'+geo_file+'----------------')
    except:
        print(output_dir)
        print('Could not glob an aims out and geometry.in - skipping')
        continue

    friction_atoms = []
    raw_tensor = []
    ptensor = False
    pvecs = False
    raw_vecs = []
    with open(output, "r") as lines:
        for i,line in enumerate(lines):

            if 'Found friction request for atom' in line:
                friction_atoms.append(int(line.split()[-1])-1)

            elif 'Printing Friction Tensor in 1/ps' in line:
                ptensor = True
                continue
            elif 'END Printing Friction Tensor' in line:
                ptensor = False

            elif ptensor:
                raw_tensor.append(line)

            elif 'Diagonalized Lifetime Eigvecs' in line:
                pvecs = True

            elif '|------------------------------------' in line or "Starting with the calculation of the polarization:" in line:
                pvecs = False

            elif 'Final output of selected total energy values:' in line:
                raw_vecs.pop()
                raw_vecs.pop()
                break
            
            elif pvecs:
                raw_vecs.append(line)

    print('Friction atoms:')
    print(friction_atoms)
    if not raw_tensor:
        print('No tensor read! - skipping')
        continue
    friction_tensor,ndim = build_tensor(raw_tensor)
    friction_vecs, ndim = build_tensor(raw_vecs)

    atoms = read(geo_file)
    modes = vjp.get_modes(atoms,friction_atoms,mode=2)

    print('---------Modes---------')
    print(modes)
    print('--------------------')

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

    np.savetxt(output_dir+'/projected_tensor.out',A,delimiter='  ',fmt='%.4e')