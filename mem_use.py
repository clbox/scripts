import numpy as np
import mem_energy_loss as mel
from ase.db import connect
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from ase.units import _hbar, J, s, fs


file_range = np.arange(1,952).tolist()
raw_data,bins= mel.Parse_memory_kernels('./calcs',file_range,read=True) #parse all memory data, choose range of datapoints parse 

#assumes in directories ./dir/n where n = 1,2,3 etc

con = connect('database.db') #connection to database that stores atoms objects and velocities

cutoffs = np.linspace(0.4,2.4,4) #all cutoff energies in eV investigating
#cutoffs = np.array([2.4])

friction_indices = [64,65] #indices of friction atoms

time_step = 2 #fs nuclear time step

mem_cutoff = 5 #fs dont include more the X fs back in time in the memory integral

pp = mel.Postprocessed_memory(bins,raw_data,cutoffs,mem_cutoff,friction_indices,time_step,con)

nm_work = pp.calculate_work()

nm_forces = pp.calculate_friction_force()

fig_path = './figures/'

dimension = len(friction_indices)*3

e=0
fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        ax[i,j].plot(pp.new_bins,(pp.new_data[0,e,:]),label='interpolated')
        ax[i,j].plot(pp.bins, pp.raw_data[0,e,:],label='raw_data')
        e+=1
    
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'freq_compare.pdf')

ts=0
e=0
fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
        ax[i,j].plot(pp.times_list[-1]/fs,(pp.eta_bar_t_list[-1])[ts,e,:],label=str(ts))
        e+=1
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_bar.pdf')

e=0
fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
        ax[i,j].plot(pp.times_up_list[-1]/fs,(pp.eta_bar_inter_list[-1])[ts,e,:],label=str(ts))
        e+=1
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_bar_inter.pdf')

times_up = pp.times_up_list[-1]
times_up = times_up[times_up >= 0.0]
times_up /= fs
e=0
fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
        ax[i,j].plot(times_up,(pp.eta_t_list[-1])[ts,e,:],label=str(ts))
        e+=1
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_t.pdf')


fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
old_time_scale = np.linspace(0,(pp.steps-1)*pp.time_step,pp.steps)
for co in range(len(cutoffs)):
    times_up = pp.times_up_list[co]
    times_up = times_up[times_up >= 0.0]
    times_up /= fs
    e=0
    for i in range(dimension):
        for j in range(i,dimension):
            eta_intgls = []
            for ts in range(len(old_time_scale)):
                eta_intgls.append(np.trapz((pp.eta_t_list[co])[ts,e,:],times_up))
            ax[i,j].plot(old_time_scale/fs,eta_intgls[:],label='CO: ' + str(cutoffs[co]))
            e+=1
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.savefig(fig_path+'eta_t_integral.pdf')


fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for co in range(len(cutoffs)):
    times_up = pp.times_up_list[co]
    times_up = times_up[times_up >= 0.0]
    times_up /= fs
    e=0
    for i in range(dimension):
        for j in range(i,dimension):
            #for ts in range(pp.steps):
            ax[i,j].plot(times_up,(pp.eta_t_list[co])[ts,e,:],label='CO: ' + str(cutoffs[co])+'TS: '+ str(ts))
            ax[i,j].set_xlim(0,15)
            e+=1
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.savefig(fig_path+'eta_t_close.pdf')

fig, ax = plt.subplots(1,1)
for atom in range(len(friction_indices)):       
    for cart in range(3):
        ax.plot(pp.inter_time_scale/fs,pp.velocities_inter[:,atom,cart],label=str(atom)+str(cart))
ax.legend()
fig.savefig(fig_path+'velocities_inter.pdf')

fig, ax = plt.subplots(1,1)
for atom in range(len(friction_indices)):       
    for cart in range(3):
        ax.plot(pp.inter_time_scale/fs,nm_forces[-1,:,atom,cart],label=str(atom)+str(cart))
ax.legend()
fig.savefig(fig_path+'nm_forces.pdf')


fig, ax = plt.subplots(1,1)
for co in range(len(cutoffs)):
    ax.plot(pp.inter_time_scale/fs,nm_work[co,:],label=str(cutoffs[co]))
ax.legend()
fig.savefig(fig_path+'nm_work.pdf')


fig, ax = plt.subplots(1,1)
for co in range(len(cutoffs)):
    ax.plot(pp.inter_time_scale/fs,np.cumsum(nm_work[co,:]),label=str(cutoffs[co]))
    #ax.plot(pp.inter_time_scale/fs,cumtrapz(nm_work[co,:],pp.inter_time_scale,initial=0),label=str(cutoffs[co]))
ax.legend()
fig.savefig(fig_path+'nm_cum_work.pdf')



#print(pp.inter_time_scale)
#print(nm_work)
datafile_path = fig_path+"nm_cumwork.txt"
cs_work = np.cumsum(nm_work[-1,:])
data = np.hstack(((pp.inter_time_scale/fs)[:,None],(cs_work[:,None])))
with open(datafile_path, 'w+') as datafile_id:
    datafile_id.write('Cutoff / eV =  {:.3f} \n'.format(cutoffs[co]))
    datafile_id.write('Time scale / fs, Work / eV \n')
    np.savetxt(datafile_id, data, fmt=['%.6f','%.6f'])
