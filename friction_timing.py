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
    with open(filename,'r') as f:

        for line in f:

            if 'Time     ' in line and read_time:
                print(line)
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


        ground_state_time = times[1]-times[0]
        finite_difference_time = times[-2]-times[1]
        tensor_time = times[-1] - times[-2]

        print(output_dir)
        print(ground_state_time)
        print(finite_difference_time)
        print(tensor_time)