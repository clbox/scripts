import numpy as np
import glob
import sys
import os
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from datetime import datetime

filenames = sys.argv[1:]

for filename in filenames:

    output_dir = os.path.dirname(filename)
    times = []
    dates = []
    read_time = True
    finite_scfs = []
    fd_scf = 0
    with open(filename,'r') as f:

        for line in f:

            if 'Begin self-consistency iteration' in line:
                ground_state_scf = int(line.split()[-1])

            if '  SCF   ' in line:

                fd_scf = int(line.split()[1])

            if ' Re-initialization.' in line and fd_scf > 0:
                finite_scfs.append(fd_scf)


            if 'Time     ' in line and read_time:
                time = line.split()[-1]
                date = (line.replace(',','').split()[-4])

                date_time = date + ' ' + time

                date_time = datetime.strptime(date_time,'%Y%m%d %H%M%S.%f')
                times.append(date_time)
                # dates.append(date)
                read_time = False

            if 'Begin self-consistency loop: Re-initialization.' in line or 'Leaving FHI-aims.' in line:
                read_time = True


        #first entry is calc start time
        #Second enetry is fric start time
        #penultimate is friction construct time
        #Last is end calc time

        print(ground_state_scf)
        print(finite_scfs)

        scfs = np.array([ground_state_scf,sum(finite_scfs),-1])

        ground_state_time = times[1]-times[0]
        finite_difference_time = times[-2]-times[1]
        tensor_time = times[-1] - times[-2]

        key_times = np.array((ground_state_time.total_seconds(),finite_difference_time.total_seconds(),tensor_time.total_seconds()))

        np.savetxt(output_dir+'/timing.txt',np.c_[key_times,scfs],fmt='%10.5f', header ='Time / s , N_cycles' )
        print(output_dir)
        print(ground_state_time)
        print(finite_difference_time)
        print(tensor_time)
