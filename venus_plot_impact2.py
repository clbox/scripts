import numpy as np
import sys
import os
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)
from matplotlib.gridspec import GridSpec

orientation = sys.argv[1] #'n' first, 'o' first or 'iso' tropic

filenames = sys.argv[2:]
for i,filename in enumerate(filenames):
    trapped = False
    print(filename)
    ntrajs = 0
    v_f=-1
    with open(filename) as f:
        if 'trapped' in filename:
            trapped = True
        else:
            v_f = int(filename.replace('.dat',''))

        for line in f:

            if 'Trajectory' in line:
                ntrajs += 1
            elif 'Initial' in line:
                #NOT trapped
                numbers = line.replace(',',' ')
                i_v = float(numbers.split()[8])
                i_r = float(numbers.split()[9])
                i_t = float(numbers.split()[10])
            elif 'Final' in line:
                #NOT trapped
                numbers = line.replace(',',' ')
                f_v = float(numbers.split()[8])
                f_r = float(numbers.split()[9])
                f_t = float(numbers.split()[10])

            elif 'Impact 0:' in line:
                pos1 = [float(a) for a in line.split()[2:]]

            elif 'Impact 1:' in line:
                pos2 = [float(a) for a in line.split()[2:]]
            
            elif 'Lifetime' in line:
                numbers = line.replace(',','')
                lifetime = float(numbers.split()[2])
                scat = float(numbers.split()[7])
                jf = int(numbers.split()[-1])




        
