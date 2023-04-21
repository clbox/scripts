import numpy as np
from ase.units import _hbar, J, s, fs
from ase import Atoms
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)
import sys
import scipy
from scipy import optimize



filenames = sys.argv[1:]


def standard_dev(x,x0,n):


    num = np.sum((x-x0)**2)

    s = np.sqrt(num/n)
    
    s = s / np.sqrt(n)
    return s





#PLOT ENERGY REDISTRIBUTION

outfile = 'raw_nfs.csv'



with open(outfile,'w') as out:

    for i,filename in enumerate(filenames):

        bin_nf = int(filename.replace('.dat',''))
        print(filename)
        with open(filename) as f:

            raw_nfs = []

            for line in f:
                if 'raw' in line:

                    nf = float(line.split()[1])
                    raw_nfs.append(nf)


                    out.write(str(nf)+',')
               
            out.write('\n')

        s = standard_dev(np.array(raw_nfs),bin_nf,len(raw_nfs))
        print(s)










