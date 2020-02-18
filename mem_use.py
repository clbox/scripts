import numpy as np
import mem_energy_loss as mel
from ase.db import connect

file_range = np.arange(1,20).tolist()
raw_data,bins,dimension = mel.Parse_memory_kernels('./calcs',file_range) #parse all memory data, choose range of datapoints parse 

#assumes in directories ./dir/n where n = 1,2,3 etc

con = connect('database.db') #connection to database that stores atoms objects and velocities

cutoffs = np.linspace(0.2,2.5,3) #all cutoff energies in eV investigating

friction_indices = [16,17,18,19] #indices of friction atoms

time_step = 2 #fs nuclear time step

mem_cutoff = 40 #fs dont include more the X fs back in time in the memory integral

pp = mel.Postprocessed_memory(bins,raw_data,cutoffs,40,friction_indices,time_step,con)

nm_work = pp.calculate_work()

nm_forces = pp.calculate_friction_force()