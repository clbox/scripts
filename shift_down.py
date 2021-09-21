from ase.io import read,write
import sys


files = sys.argv[1:]

for file in files:
    print(file)

    slab = read(file)

    z = slab.positions[0,2]
    slab.positions[:,2] -= z

    write(file,slab)

    slab = None