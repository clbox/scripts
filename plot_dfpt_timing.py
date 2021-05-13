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
colours = ['black','red','blue','orange','brown','darkgreen','purple']
marker_size = 3

for i,filename in enumerate(filenames):
    x_vals = []
    y_vals = []
    s1_vals = []
    h1_vals = []
    u1_vals = []
    rho1_vals = []
    v1_vals = []
    dm1_vals = []
    use_dfpt = False
    with open(filename,'r') as f:
        for line in f:
            if 'Filename' in line:
                current_x_val = None
                current_y_val = None


             #Assuming all same in one file
            if 'Scalapack' in line:
                if 'False' in line:
                    label = 'LAPACK'
                if 'True' in line:
                    label = 'ScaLAPACK'

            if x_string in line:
                current_x_val = int(line.split()[-1])
            if y_string in line:
                current_y_val = float(line.split()[-1])

            if 'Total_time_ ' in line:
                x_vals.append(current_x_val)
                y_vals.append(current_y_val)

            if 'first_order_DM' in line:
                dm1_vals.append(float(line.split()[-1]))

            if 'Total_time_first_order_density' in line:
                rho1_vals.append(float(line.split()[-1]))

            if 'first_order_potential' in line:
                v1_vals.append(float(line.split()[-1]))

            if 'first_order_H' in line:
                h1_vals.append(float(line.split()[-1]))

            if 'Sternheimer' in line:
                u1_vals.append(float(line.split()[-1]))

            if 'first_order_S' in line:
                calc_s1 = True
                s1_vals.append(float(line.split()[-1]))
            

    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)

    index  = np.argsort(x_vals)
    x_vals = x_vals[index]
    y_vals = y_vals[index]



    ax.plot(x_vals,y_vals,
        mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
        label='Total time',marker=markers[i],color=colours[0],linestyle=linestyles[i])


    ax.plot(x_vals,h1_vals,mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
        label='H(1)',marker=markers[i],color=colours[1],linestyle=linestyles[i])

    ax.plot(x_vals,u1_vals,mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
        label='U(1)',marker=markers[i],color=colours[2],linestyle=linestyles[i])

    ax.plot(x_vals,v1_vals,mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
        label='V(1)',marker=markers[i],color=colours[3],linestyle=linestyles[i])

    ax.plot(x_vals,rho1_vals,mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
        label=r'$\rho$(1)',marker=markers[i],color=colours[4],linestyle=linestyles[i])

    ax.plot(x_vals,dm1_vals,mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
        label='DM(1)',marker=markers[i],color=colours[5],linestyle=linestyles[i])

    if calc_s1:
        ax.plot(x_vals,s1_vals,mfc='none',linewidth=0.4,markersize=marker_size,markeredgewidth=0.4,
            label='S(1)',marker=markers[i],color=colours[6],linestyle=linestyles[i])




ax.set_xlim(left=0)
ax.legend()

ax.set_yscale('log')

#ax.set_xlabel('Number of atoms')
ax.set_xlabel('Number of processes')
ax.set_ylabel('Total time / s')

fig.set_figheight(5.)
fig.set_figwidth(5.)
fig.savefig('timing.pdf',transparent=False,bbox_inches='tight')
