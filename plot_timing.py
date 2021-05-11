import sys
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


#Reads output from timing.py

filename = sys.argv[1]

x_string = 'N_atoms'
#x_string = 'N_tasks'

y_string = 'Total_time'

x_vals = []
y_vals = []

with open(filename,'r') as f:
    for line in f:
        if 'Filename' in line:
            current_x_val = None
            current_y_val = None

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

fig, ax = plt.subplots(1, 1)

ax.plot(x_vals,y_vals,
    linestyle='-',marker='o',color='black',mfc='none',linewidth=0.4,markersize=5,markeredgewidth=0.4)

ax.set_yscale('log')

ax.set_xlabel('Number of atoms')
#ax.set_xlabel('Number of processes')
ax.set_ylabel('Total time / s')

fig.set_figheight(5.)
fig.set_figwidth(5.)
fig.savefig('timing.pdf',transparent=True,bbox_inches='tight')
