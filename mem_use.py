import numpy as np
import mem_energy_loss as mel
from ase.db import connect
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from ase.units import _hbar, J, s, fs
from scipy.integrate import cumtrapz


file_range = np.arange(1,101).tolist()
raw_data,bins,dimension = mel.Parse_memory_kernels('./calcs',file_range)#,read=True) #parse all memory data, choose range of datapoints parse 

#assumes in directories ./dir/n where n = 1,2,3 etc

con = connect('database.db') #connection to database that stores atoms objects and velocities

#cutoffs = np.linspace(0.2,2.4,2) #all cutoff energies in eV investigating
cutoffs = [2.4]

friction_indices = [64,65] #indices of friction atoms

time_step = 2 #fs nuclear time step

mem_cutoff = 40 #fs dont include more the X fs back in time in the memory integral

pp = mel.Postprocessed_memory(bins,raw_data,cutoffs,mem_cutoff,friction_indices,time_step,con)

nm_work = pp.calculate_work()

nm_forces = pp.calculate_friction_force()

fig_path = './figures/'

fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        ax[i,j].plot(pp.new_bins,(pp.new_data[0,i+j,:]),label='interpolated')
        ax[i,j].plot(pp.bins, pp.raw_data[0,i+j,:],label='raw_data')
    
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'freq_compare.pdf')

ts=0
fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
        ax[i,j].plot(pp.times_list[-1]/fs,(pp.eta_bar_t_list[-1])[ts,i+j,:],label=str(ts))
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_bar.pdf')


fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
            ax[i,j].plot(pp.times_up_list[-1]/fs,(pp.eta_bar_inter_list[-1])[ts,i+j,:],label=str(ts))
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_bar_inter.pdf')


fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
        ax[i,j].plot(pp.times_up_list[-1]/fs,(pp.eta_bar_inter_list[-1])[ts,i+j,:],label=str(ts))
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_inter.pdf')

fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for i in range(dimension):
    for j in range(i,dimension):
        #for ts in range(pp.steps):
        ax[i,j].plot(pp.times_up_list[-1]/fs,(pp.eta_t_list[-1])[ts,i+j,:],label=str(ts))
ax[0,0].legend()
fig.set_figheight(20)
fig.set_figwidth(20)
fig.savefig(fig_path+'eta_t.pdf')

fig, ax = plt.subplots(dimension,dimension,sharex='all', sharey='all')
for co in range(len(cutoffs)):
    for i in range(dimension):
        for j in range(i,dimension):
            #for ts in range(pp.steps):
            ax[i,j].plot(pp.times_up_list[co]/fs,(pp.eta_t_list[co])[ts,i+j,:],label='CO: ' + str(cutoffs[co])+'TS: '+ str(ts))
            ax[i,j].set_xlim(0,5)
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
    ax.plot(pp.inter_time_scale/fs,cumtrapz(nm_work[co,:],pp.inter_time_scale/fs,initial=0),label=str(cutoffs[co]))
ax.legend()
fig.savefig(fig_path+'nm_cum_work.pdf')

#print(pp.inter_time_scale)
#print(nm_work)

np.savetxt('nm_work.txt',(pp.inter_time_scale,nm_work))