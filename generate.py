from ase.build import molecule
from ase.io import read,write
from ase.build import surface, add_adsorbate, fcc111, hcp0001

# Build slab
slab = hcp0001('Ru',a=2.658,size=(2,2,6),vacuum=50)

# Add adsorbate
atoms = molecule('CO')
h = 3 #height ,angstrom
add_adsorbate(slab, atoms, h, 'ontop') #mol_index to pick which atom is ontop

# Write file
write('geometry.in', slab, format='aims')
