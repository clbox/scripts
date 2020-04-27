from ase import Atoms
from ase.io import read,write
from ase.io.aims import write_aims
import numpy as np


system = read('geometry.in')

default_positions  = system.get_positions()

surface_height = np.max(system.get_positions()[2:,2])


for i,height in enumerate(np.linspace(0.1,7.1,15)):
    print(height)
    new_positions = default_positions
    new_positions[0,2] = height + surface_height
    new_positions[1,2] = height + surface_height
    system.set_positions(new_positions)

    write_aims('geometry_{:02d}.in'.format(i),system,friction_atoms=[0,1])


