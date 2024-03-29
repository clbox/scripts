import sys
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


#Reads output from timing.py

filenames = sys.argv[1:]

fig, ax = plt.subplots(1, 1)

#x_string = 'N_atoms'
x_string = 'N_tasks'

y_string = 'Total_time'

markers = ['o','^','s','.']
linestyles = ['-',':','-.','--']

for i,filename in enumerate(filenames):
    x_vals = []
    y_vals = []
    use_dfpt = False
    marker = 'o'
    with open(filename,'r') as f:
        for line in f:
            if 'Filename' in line:
                current_x_val = None
                current_y_val = None

            if 'DFPT' in line:
                use_dfpt = True

             #Assuming all same in one file
            if 'Scalapack' in line:
                if 'False' in line:
                    label = 'LAPACK'
                if 'True' in line:
                    label = 'ScaLAPACK'
                    marker = '^'

            if x_string in line:
                current_x_val = int(line.split()[-1])
            if y_string in line:
                current_y_val = float(line.split()[-1])

            if 'Total_time_ ' in line:
                x_vals.append(current_x_val)
                y_vals.append(current_y_val)

    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)

    index  = np.argsort(x_vals)
    x_vals = x_vals[index]
    y_vals = y_vals[index]


    if use_dfpt:
        label = 'DFPT '+label
        color = 'blue'
    else:
        label = 'FD '+label
        color = 'black'

    if 'really_tight' in filename:
        label = label + ' Really_tight'
        linestyle = '--'
    else:
        linestyle = '-'
        
    ax.plot(x_vals,y_vals,
        mfc='none',linewidth=0.4,markersize=3,markeredgewidth=0.4,
        label=label,marker=marker,color=color,linestyle=linestyle)

ax.set_xlim(left=0)
ax.set_ylim(bottom=1e3,top=1e5)
ax.legend()

ax.set_yscale('log')

#ax.set_xlabel('Number of atoms')
ax.set_xlabel('Number of processes')
ax.set_ylabel('Total time / s')

fig.set_figheight(5.)
fig.set_figwidth(5.)
fig.savefig('timing.pdf',transparent=False,bbox_inches='tight')
