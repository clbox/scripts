import numpy as np
from ase.units import _hbar, J, s, fs, kcal, mol
from ase import Atoms
import glob
import sys
import os
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


filenames = sys.argv[1:]
outfile = 'energies.out'
for i, filename in enumerate(filenames):

    with open(filename) as f:
        for line in f:
            if '************************FRICTION**********************************' in line:
                break
            if '   | Total energy corrected       ' in line:
                e_corr = float(line.split()[5])
                break

    with open(outfile,'a') as f:
        f.write(str(i)+'   '+e_corr)
    

