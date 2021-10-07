import numpy as np
import sys
import os
import glob
from pathlib import Path
import matplotlib
from numpy.testing._private.utils import jiffies
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from ase.io import read, write


# User gives geometry file
input_file = sys.argv[1]

# User gives atoms for structure
surface_element = 'Cu'
adsorbate_elements = ['C','O']

# Read geometry to atoms object
structure = read(input_file)

# Find average height of top layer of surface
surface_indices = [i for i, e in enumerate(structure) if e.symbol == surface_element]
surface_atoms = structure[surface_indices]
max_surface_atom_height = np.max(surface_atoms.positions[:,2])
tol = 0.3
top_layer_indices = np.argwhere(surface_atoms.positions[:,2]>max_surface_atom_height-tol)
average_top_layer_z = np.mean(surface_atoms.positions[top_layer_indices,2])
print('Average top layer z / Ang : ' + str(average_top_layer_z))

# Print all heights of each adsorbate atom
average_adsorbate_heights = []
for element in adsorbate_elements:
    element_indices = [i for i, e in enumerate(structure) if e.symbol == element]
    element_atoms = structure[element_indices]
    average_element_z = np.mean(element_atoms.positions[:,2])

    print('Printing z coordinates for element: ' + str(element) + ' / Angstrom')
    print(element_atoms.positions[:,2])

    print('Average ' +str(element) + ' z / Ang : '+str(average_element_z))

    # Print average height of each adsorbate atom
    average_adsorbate_height=average_element_z-average_top_layer_z
    print('Average ' +str(element) + ' height / Ang : '+str(average_element_z-average_top_layer_z))
    average_adsorbate_heights.append(average_adsorbate_height)


for i,element1 in enumerate(adsorbate_elements):
    for j,element2 in enumerate(adsorbate_elements):

        if element1==element2:
            continue

        if j < i:
            continue

        print('Average height difference between ' + str(element1) + ' and ' + str(element2) +' is ' + str(np.abs(average_adsorbate_heights[j]-average_adsorbate_heights[i])))
