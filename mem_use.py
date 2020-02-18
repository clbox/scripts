import mem_energy_loss as me


raw_data,bins,dimension = me.Parse_memory_kernels('./calcs',np.arange(1,20)) #parse all memory data, choose range of datapoints parse 
#assumes in directories ./dir/n where n = 1,2,3 etc

con = connect('database.db') #connection to database that stores atoms objects and velocities

cutoffs = np.linspace(0.2,2.5,3) #all cutoff energies in eV investigating

friction_indices = [16,17,18,19] #indices of friction atoms

time_step = 2 #fs nuclear time step


pp = me.Postprocessed_memory(bins,raw_data,cutoffs,40,friction_indices,time_step,con)

nm_work = pp.calculate_work()

nm_forces = pp.calculate_friction_force()