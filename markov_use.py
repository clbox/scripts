import numpy as np
import mem_energy_loss as mel
from ase.db import connect
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from ase.units import _hbar, J, s, fs
ps = fs*1000
import os


con = connect('database.db') #connection to database that stores atoms objects and velocities

friction_indices = [64,65] #indices of friction atoms

time_step = 2 #fs nuclear time step

steps = 952

key = 'ft1' 

pp = mel.Postprocessed_markov(steps,friction_indices,time_step,con,key)

m_work = pp.calculate_work()

m_forces = pp.calculate_friction_force()

fig_path = './markov_figures/'

dimension = len(friction_indices)*3

all_velocities = pp.all_velocities

fig, ax = plt.subplots(1,1)
for atom in range(len(friction_indices)):       
    for cart in range(3):
        ax.plot(pp.time_scale/fs,all_velocities[:,atom,cart]*fs,label=str(atom)+str(cart))
ax.legend()
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'Velocity /  $\AA$ fs$^{-1}$', va='center', rotation='vertical',fontsize=15)
fig.savefig(fig_path+'all_velocities.pdf')

all_tensors = pp.get_all_tensors()
fig, ax = plt.subplots(1,1)
labels = ['C_x','C_y','C_z','O_x','O_y','O_z']
for i in range(dimension):
    for j in range(i,dimension):
        if j != i:
            continue
        ax.plot(pp.time_scale/fs,all_tensors[:,i,j]*ps,label=str(labels[i])+' '+str(labels[j]))
ax.legend()
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'$\Lambda_{ij}(0) $/ ps$^{-1}$', va='center', rotation='vertical',fontsize=15)
fig.savefig(fig_path+'all_tensors.pdf')


fig, ax = plt.subplots(1,1)
for atom in range(len(friction_indices)):       
    for cart in range(3):
        ax.plot(pp.time_scale/fs,m_forces[:,atom,cart],label=str(atom)+str(cart))
ax.legend()
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'Friction force / eV $\AA^{-1}$', va='center', rotation='vertical',fontsize=15)
fig.savefig(fig_path+'m_forces.pdf')


fig, ax = plt.subplots(1,1)
ax.plot(pp.time_scale/fs,m_work[:])
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Work / eV', va='center', rotation='vertical',fontsize=15)
fig.savefig(fig_path+'m_work.pdf')


fig, ax = plt.subplots(1,1)
ax.plot(pp.time_scale/fs,np.cumsum(m_work[:]))
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Cumulative work / eV', va='center', rotation='vertical',fontsize=15)
fig.savefig(fig_path+'m_cum_work.pdf')



datafile_path = fig_path+"m_cum_work.txt"
try:
    os.remove(datafile_path)
except OSError:
    pass
cs_work = np.cumsum(m_work[:])
data = np.hstack(((pp.time_scale/fs)[:,None],(cs_work[:,None])))
with open(datafile_path, 'w+') as datafile_id:
    datafile_id.write('Time scale / fs, Work / eV \n')
    np.savetxt(datafile_id, data, fmt=['%.6f','%.6f'])