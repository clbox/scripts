from ase.build import molecule
from ase.io import read,write
from ase.build import surface, add_adsorbate, fcc111, hcp0001, sort
import numpy as np
import os

# Define slab
substrate = 'Cu'
min_layers = 3
# Parse structure
full_structure = read('geometry.in')

# Sort structure depending on z height
idx = np.argsort(full_structure.positions[:,2])
print(full_structure.positions[:,2])
print(idx)
full_structure = full_structure[idx]
#full_structure = sort(full_structure,tags=idx)
print(full_structure.positions[:,2])

# Find atoms per layer (assume same for now)
bottom_layer_height = full_structure.positions[0,2]
atoms_per_layer = np.count_nonzero(full_structure.positions[:,2]==bottom_layer_height)
print('I think there is ' + str(atoms_per_layer) + ' atoms per layer')

# Find number of layers
substrate_index = [i for i, e in enumerate(full_structure) if e.symbol == substrate]
A = full_structure.positions[substrate_index,2]
tol = 0.3
number_of_layers = len(A[~(np.triu(np.abs(A[:,None] - A) <= tol,1)).any(0)])
print('I think there is ' + str(number_of_layers) + ' layers')

# Interlayer distance
inter_layer_distance = full_structure.positions[atoms_per_layer,2]

# Remove layer , reduce cell height, make subdirectory, write file to subdirectory
for l in range(number_of_layers-min_layers):

    # Remove layer of atoms
    for i in range(atoms_per_layer):
        full_structure.pop(0)


    # Create subdirectory
    sub_dir = '{:02d}'.format(number_of_layers-(l+1))
    os.mkdir(sub_dir)

    # Shift down
    z = full_structure.positions[0,2]
    full_structure.positions[:,2] -= z

    # Shrink cell height accordingly
    cell = full_structure.get_cell()
    cell[2,2] -= inter_layer_distance
    full_structure.set_cell(cell)

    # Write file
    write(sub_dir+'/geometry.in', full_structure, format='aims')

    

