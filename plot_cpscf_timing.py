import sys
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

matplotlib.rcParams['font.sans-serif'] = "Arial"
#Reads output from timing.py

filenames = sys.argv[1:]

fig, ax = plt.subplots(1, 1)

#x_string = 'N_atoms'
x_string = 'N_tasks'
y_string = 'Total_time'

markers = ['o','.','^','s']
linestyles = ['-',':','-.','--']
colours = ['red','blue','orange','brown','darkgreen','purple']
marker_size = 3




labels = ['CPSCF legacy','CPSCF (v0.01)']

for m,mode in enumerate(['legacy','v001']):
    n_atoms_per_file = []
    avg_iter_time_per_file = []
    avg_scf_iter_time_per_file = []
    for i,filename in enumerate(filenames):

        if mode not in filename:
            continue

        iteration_times = []
        scf_iteration_times = []
        read_timing = False
        read_scf_timing = False

        with open(filename,'r') as f:
            for line in f:

                if '| Number of atoms' in line:
                    n_atoms = int(line.split()[-1])
                    n_atoms_per_file.append(n_atoms)
                    continue

                if 'End CPSCF iteration # ' in line:
                    read_timing = True
                    continue

                if 'End self-consistency iteration #' in line:
                    read_scf_timing = True
                    continue

                if '-------------------------------------' in line:
                    read_timing = False
                    read_scf_timing = False
                    continue

                if '| Time for this iteration' in line and read_scf_timing:
                    time_for_iteration = float(line.split()[-4])
                    scf_iteration_times.append(time_for_iteration)
                    continue

                if '| Time for this iteration' in line and read_timing:
                    time_for_iteration = float(line.split()[-4])
                    iteration_times.append(time_for_iteration)
                    continue

        iteration_times = np.array(iteration_times)
        scf_iteration_times = np.array(scf_iteration_times)

        average_iter_time = np.mean(iteration_times)
        average_scf_iter_time = np.mean(scf_iteration_times)

        avg_iter_time_per_file.append(average_iter_time)
        avg_scf_iter_time_per_file.append(average_scf_iter_time)


        
    

    n_atoms_per_file = np.array(n_atoms_per_file)
    avg_iter_time_per_file = np.array(avg_iter_time_per_file)
    avg_scf_iter_time_per_file = np.array(avg_scf_iter_time_per_file)

    idx = np.argsort(n_atoms_per_file)

    n_atoms_per_file = n_atoms_per_file[idx]
    avg_iter_time_per_file = avg_iter_time_per_file[idx]
    avg_scf_iter_time_per_file = avg_scf_iter_time_per_file[idx]


    ax.plot(n_atoms_per_file,avg_iter_time_per_file,linestyle='-',linewidth=0.7,marker=markers[m],mfc='none',color=colours[m],label=labels[m])

    if 'legacy' in mode:
        continue

    ax.plot(n_atoms_per_file,avg_scf_iter_time_per_file,linestyle='-',linewidth=0.7,marker='^',mfc='none',color='black',label='SCF')



# x = np.linspace(0,np.max(n_atoms_per_file),200)
# y = x**3
# ax.plot(x,y,linestyle='--',linewidth=0.1)
ax.set_xscale('log')
# ax.set_xticks(n_atoms_per_file)
# ax.set_xticklabels(n_atoms_per_file)

ax.set_xticks([10,50,100,250,500])
ax.set_xticklabels([10,50,100,250,500])

ax.set_xlabel('Number of atoms')
ax.set_ylabel('Average time per iteration / s')
ax.set_yscale('log')

ax.legend(fancybox=True,framealpha=0)

fig.set_figheight(3.5)
fig.set_figwidth(4.)
fig.savefig('timing.pdf',transparent=False,bbox_inches='tight')



            