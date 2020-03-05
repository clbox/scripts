import numpy as np
import mem_energy_loss as mel
from ase.db import connect
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from ase.units import _hbar, J, s, fs
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
        ax.plot(pp.time_scale/fs,all_velocities[:,atom,cart],label=str(atom)+str(cart))
ax.legend()
fig.savefig(fig_path+'all_velocities.pdf')

fig, ax = plt.subplots(1,1)
for atom in range(len(friction_indices)):       
    for cart in range(3):
        ax.plot(pp.time_scale/fs,m_forces[:,atom,cart],label=str(atom)+str(cart))
ax.legend()
fig.savefig(fig_path+'m_forces.pdf')


fig, ax = plt.subplots(1,1)
ax.plot(pp.time_scale/fs,m_work[:])
fig.savefig(fig_path+'m_work.pdf')


fig, ax = plt.subplots(1,1)
ax.plot(pp.time_scale/fs,np.cumsum(m_work[:]))
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