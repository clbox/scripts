from ase.build import molecule
from ase.io import read,write
from ase.build import surface, add_adsorbate, fcc111, hcp0001
import numpy as np
import os
import sys

# Read geometry.in
input_file_name = sys.argv[1]
structure = read(input_file_name)


# Set QE filename
output_file_name=input_file_name.replace('geometry','scf')

# Reduce vacuum as costly for PW
cell = structure.get_cell()
cell[2,2] -= 120.0
structure.set_cell(cell)

# Write QE file
print('Written to: ', output_file_name)
write(output_file_name, structure, format='espresso-in')
