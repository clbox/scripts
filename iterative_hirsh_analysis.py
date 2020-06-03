import numpy as np
import sys
import os
#Reads in several aims.out paths, outputs final hirsh charges to txt file in that path as well as sum of charges


paths = sys.argv[1:]


for path in paths:
    output_dir = os.path.dirname(path)


    with open(path,'r') as f:
        for line in f:
            if '| Number of atoms' in line:
                natoms = int(line.split()[-1])
            
            if 'Hirshfeld-I iterations' in line:
                final_charges = []

            if '|   Hirshfeld charge ' in line:
                final_charges.append(float(line.split()[-1]))



    final_charges = np.array(final_charges)


    np.savetxt(output_dir+'/hirsh.txt',final_charges,delimiter='  ',fmt='%.4e')
    print(np.sum(final_charges))