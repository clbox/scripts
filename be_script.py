from ase import Atoms
from ase.io import read,write
from ase.io.aims import write_aims
import numpy as np


system = read('geometry.in')

atom1 = 0
atom2 = 1

default_positions  = system.get_positions()

masses = system.get_masses()
m1 = masses[atom1]
m2 = masses[atom2]
com = (m1*default_positions[atom1] + m2*default_positions[atom2])/(m1+m2)
a1_relative_height = default_positions[atom1,2] - com[2]
a2_relative_height = default_positions[atom2,2] - com[2]

surface_height = np.max(system.get_positions()[2:,2])

#relative_height = default_positions[atom1,2] -  default_positions[atom2,2]


for i,height in enumerate(np.linspace(0.1,7.1,15)):
    print(height)
    new_positions = default_positions
    new_positions[atom1,2] = height + surface_height + a1_relative_height
    new_positions[atom2,2] = height + surface_height + a2_relative_height
    system.set_positions(new_positions)

    write_aims('geometry_{:02d}.in'.format(i),system,friction_atoms=[0,1])


