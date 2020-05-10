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
import venus_traj_project as vjp 
from ase.visualize import view



#read in out26, if trajectory no not in there then append to list of trapped

#get atoms_list for each trapped

#and plot simple graph showing o and n positions?

#TODO: plot the impact geometry . I.e when COM is lowest, plot geometry of NO (x vs z?)
#TODO: add support for choosing final vib state, not just trapped
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')


label = True

mode = sys.argv[1]

impact_geo = []

traj_nos = sys.argv[2:]
for traj_no in traj_nos:
    O_pos = []
    N_pos = []
    traj_no = int(traj_no)
    va = vjp.venus_analysis(traj_no,[0,1],mode=2)

    if mode == 'trapped':
        trapped = va.is_trajectory_trapped()

        if not trapped:
            continue

    else:
        try:
            lifetime,Nf,Jf,scat_angle = va.parse_traj_summary()
        except:
            #trappped
            continue
        Nf = va.bin_quantum(Nf)
        if int(mode) != Nf:
            continue 


    trajectory = va.build_atoms_list()

    for step in trajectory:
        O_pos.append(step.get_positions()[0])
        N_pos.append(step.get_positions()[1])

    O_pos = np.array(O_pos)
    N_pos = np.array(N_pos)

    print(len(O_pos))
    
    ON_z_pos = np.column_stack((O_pos[0:300,2],N_pos[0:300,2]))

    idx = np.unravel_index(np.argmin(ON_z_pos, axis=None), ON_z_pos.shape)
    idx = idx[0]

    print(idx)

    print('---')

    impact_geo.append([O_pos[idx],N_pos[idx]])
    
    
    #looking for where an atom gets closest to surface



    times = np.arange(1,len(O_pos)+1,1)
    times = times #NIP = 10

    if label:
        ax.plot(times,O_pos[:,2],color='red',linewidth=0.5,label='Oxygen')
        ax.plot(times,N_pos[:,2],color='blue',linewidth=0.5,label='Nitrogen')
        label = False
    else:
        ax.plot(times,O_pos[:,2],color='red',linewidth=0.5)
        ax.plot(times,N_pos[:,2],color='blue',linewidth=0.5)



ax.set_xlim(0,600)
ax.set_ylim(0,7)
ax.legend(loc=1)

fig.set_figheight(4)
fig.set_figwidth(5)
fig.text(0.5, 0.00, r"Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'$z$ / $\AA{}$', va='center', rotation='vertical',fontsize=15)
fig.savefig('traj_plot_'+mode+'.pdf',transparent=True,bbox_inches='tight')



fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')

impact_geo = np.array(impact_geo)

ax.scatter(impact_geo[:,0,0],impact_geo[:,0,2],c='red',s=800,alpha=0.7,label='Oxygen')
ax.scatter(impact_geo[:,1,0],impact_geo[:,1,2],c='blue',s=800,alpha=0.7,label='Nitrogen')

ax.set_xlim(-5,6)
ax.set_ylim(0,7)
lgnd = ax.legend(loc=1)
lgnd.legendHandles[0]._sizes = [200]
lgnd.legendHandles[1]._sizes = [200]

fig.set_figheight(4)
fig.set_figwidth(5)
fig.text(0.5, 0.00, r"$x$ / $\AA{}$", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'$z$ / $\AA{}$', va='center', rotation='vertical',fontsize=15)
fig.savefig('impact_'+mode+'.pdf',transparent=True,bbox_inches='tight')


with open('impact_'+mode+'.txt', 'w') as outfile:
    # I'm writing a header here just for the sake of readability
    # Any line starting with "#" will be ignored by numpy.loadtxt
    outfile.write('# Array shape: {0}\n'.format(impact_geo.shape))

    # Iterating through a ndimensional array produces slices along
    # the last axis. This is equivalent to data[i,:,:] in this case
    for data_slice in impact_geo:

        # The formatting string indicates that I'm writing out
        # the values in left-justified columns 7 characters in width
        # with 2 decimal places.  
        np.savetxt(outfile, data_slice, fmt='%-7.2f')

        # Writing out a break to indicate different slices...
        outfile.write('# New traj\n')
#ax.plot(impact_geo)