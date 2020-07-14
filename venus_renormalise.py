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
                            

#Reads in N amount of states files
# Determines what vib state looking at from file directory
#Renormalises data according to experiment
#Saves file with e appended to filename


filenames = sys.argv[1:]


for i,filename in enumerate(filenames):
    exp_states = []
    exp_dist = []
    dist = np.loadtxt(filename)

    if 'v02' in os.path.abspath(filename):
        exp_range = [1,3]

    elif 'v03' in os.path.abspath(filename):
        exp_range = [1,3]

    elif 'v11' in os.path.abspath(filename):
        exp_range = [2,11]

    elif 'v16' in os.path.abspath(filename):
        exp_range = [2,16]

    else:
        print('Could not find vibrational state - continuing')
        continue

    for i,state in enumerate(dist[:,0]):
        if state < exp_range[0] or state > exp_range[1]:
            continue

        else:
            exp_states.append(state)
            exp_dist.append(dist[i,1])

    new_raw_dist = np.column_stack((exp_states, exp_dist))

    new_dist = []
    for i in range(len(new_raw_dist[:,0])):
        normed = new_raw_dist[i,1] / np.sum(new_raw_dist[:,1])
        new_dist.append(normed)



    np.savetxt(filename.replace('.txt','_e.txt'), np.c_[exp_states,new_dist],fmt='%1.3f')



