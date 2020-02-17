#Creates a sequence of geometry files from just one, with a range of bond lengths for a single selected bond
#There is probably an easier way to do this

import numpy as np
from ase.io import read, write

##################################################################
##################################################################

index_atom1 = 64  #indices of two atoms whose bond you are changing (indices in Atoms object)
index_atom2 = 65
n_atoms = 66 # number of atoms in system
drs = np.linspace(-0.15,0.15,25)  #range of bond length differences you want -ve = shrink +ve = expand

##################################################################
##################################################################

for ii,dr in enumerate(drs):
    t = 0.5*(dr**2)*(1/np.sqrt(3))
    v = np.zeros((n_atoms,3))
    if dr > 0:
        v[index_atom1,:] = [t,t,t]
        v[index_atom2,:] = [-t,-t,-t]
    elif dr < 0:
        v[index_atom1,:] = [-t,-t,-t]
        v[index_atom2,:] = [t,t,t]
    elif dr == 0:
        print('same geometry')
    ab = read('start.in')
    ab.translate(v)
    write('geometry_'+'{:02d}'.format(ii)+'.in', ab, format='aims')
