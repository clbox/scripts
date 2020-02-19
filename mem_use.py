import numpy as np
import mem_energy_loss as mel
from ase.db import connect
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from ase.units import _hbar, J, s, fs


file_range = np.arange(1,20).tolist()
raw_data,bins,dimension = mel.Parse_memory_kernels('./calcs',file_range,read=True) #parse all memory data, choose range of datapoints parse 

#assumes in directories ./dir/n where n = 1,2,3 etc

con = connect('database.db') #connection to database that stores atoms objects and velocities

cutoffs = np.linspace(0.2,2.5,3) #all cutoff energies in eV investigating

friction_indices = [16,17,18,19] #indices of friction atoms

time_step = 2 #fs nuclear time step

mem_cutoff = 40 #fs dont include more the X fs back in time in the memory integral

pp = mel.Postprocessed_memory(bins,raw_data,cutoffs,40,friction_indices,time_step,con)

nm_work = pp.calculate_work()

nm_forces = pp.calculate_friction_force()

fig, ax = plt.subplots(1,1)
ax.plot(pp.times_up_list[-1]/fs,(pp.eta_bar_inter_list[-1])[0,0,0,:])
fig.savefig('eta_inter.pdf')

fig, ax = plt.subplots(1,1)
for atom in range(len(friction_indices)):       
    for cart in range(3):
        ax.plot(pp.inter_time_scale,pp.velocities_inter[:,atom,cart],label=str(atom)+str(cart))
ax.legend()
fig.savefig('velocities_inter.pdf')





print(pp.inter_time_scale)
print(nm_work)

np.savetxt('nm_work.txt',(pp.inter_time_scale,nm_work))