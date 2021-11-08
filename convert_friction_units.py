import numpy as np
import sys
import os
from ase import Atoms

from ase.units import  J, fs 
ps = fs*1000
meV = 1000

def parse_fiction_masses(aims_file):
    friction_atoms = []
    with open(aims_file, "r") as af:
        read = False
        for line in af:
            if 'The contents of geometry.in will be repeated verbatim below' in line:
                read=True
            if 'calculate_friction .true.' in line:
                friction_atoms.append(element)
            if read:
                try:
                    element = line.split()[-1]
                except:
                    continue
    
    a = Atoms(symbols=friction_atoms)
    friction_masses = a.get_masses()

    return(friction_masses)

aims_file = sys.argv[1] #Just for atoms
friction_tensor_file = sys.argv[2]  #friction_tensor.out (mass_weighted)

#could overwrite tensor if not in .out file
output_file = friction_tensor_file.replace('.out','proper.out')

friction_masses = parse_fiction_masses(aims_file)
print(friction_masses)


mass_weighted_ft = np.loadtxt(friction_tensor_file)
ft = np.zeros_like(mass_weighted_ft)


dimension = np.shape(mass_weighted_ft)[0]
for i in range(dimension):
    i_atom = i // 3
    for j in range(dimension):
        j_atom = j // 3
        ft[i,j] = mass_weighted_ft[i,j]*np.sqrt(friction_masses[i_atom]*friction_masses[j_atom])

# ft is now not mass weighted but in units of amu ps-1
ft = ft / ps 
ft = ft * meV 
ft = ft / ps

print(ft) #meV ps Ang-2

np.savetxt(output_file,ft)