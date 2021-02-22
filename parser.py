from ase import Atoms,Atom
from ase.io import read,write
import os
import numpy as np
from ase.io.trajectory import TrajectoryWriter
import sys
from ase.io.trajectory import TrajectoryReader
from ase.db import connect


#atoms_traj = TrajectoryReader('tdpt_1.traj')
atoms_traj = read('tdpt_1_rightvel.traj@:')


lattice_vectors=[      8.85439968,      0.00000000,      0.00000000,
     -4.42719984,      7.66813500,      0.00000000,
      0.00000000,      0.00000000,     72.22960090,]

lattice_vectors = np.reshape(lattice_vectors,(3,3))

print(lattice_vectors)

W = TrajectoryWriter('tdpt_1_cell.traj')

with connect('traj_database.db') as con:
    #for i in range(len(atoms_traj)):
    for i in range(len(atoms_traj)):
        #ab.set_cell(lv)       
        

        ab = atoms_traj[i]

        
        ab.set_cell(lattice_vectors)
        
        ab.set_pbc((True,True,True))

    
        #l=ab.get_cell_lengths_and_angles()
        #convert fractional to cartesian
       # for ii in range(np.shape(atoms[i,:,:])[0]):
       #     p = atoms[i,ii,:]
       #     if ii==0:         
       #         ab.append(Atom('O',(p[0]*l[0],p[1]*l[1],p[2]*l[2])))
       #     elif ii==1:         
       #         ab.append(Atom('N',(p[0]*l[0],p[1]*l[1],p[2]*l[2])))
       #     else:                                     
       #         ab.append(Atom('Au',(p[0]*l[0],p[1]*l[1],p[2]*l[2])))
    
        #increase z lattice length
        #l[2] += 20
        #ab.set_pbc((True,True,True))
        #ab.set_cell(l)
        W.write(ab)
   # OVERWRITES DATABASE CAREFUL!!!!!
        
        con.write(ab,calculation_status=0)
