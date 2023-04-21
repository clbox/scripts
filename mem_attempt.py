#!/usr/bin/env python

import numpy as np
import scipy.sparse as sp
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                       AutoMinorLocator, MaxNLocator)
import elsi_mem_parsing as emp
import time
from ase.units import _hbar, J, s, fs, Bohr

hbar = _hbar * J * s 
ps = fs*1000

plt.style.use('clb_publication')

# one_over_hbar = 1/(6.582119569e-16*1e15) #ev fs
one_over_hbar = 1/hbar

fig, ax = plt.subplots(1, 1)


# k_weights = [0.01234568] + [0.02469136] * 40
k_weights = [0.25]*4

# Define limits of integration
max_energy_from_fermi = 1.
n_k_points = 4 #41

i_spin = 0

mem_cutoff = 10 #fs


# steps = [109,110,111] # base 0
steps = list(range(10))

time_step_delta = (steps[1]-steps[0])*0.1

time_steps = np.array(steps) * 0.1

friction_atoms = [0,1]
n_cart = 3
n_friction_atoms = len(friction_atoms)

friction_kernel = np.zeros((n_friction_atoms,n_cart,n_friction_atoms,n_cart,len(steps),len(steps)))

for i_atom in friction_atoms:
    for i_cart in range(n_cart):
        

        for i,step in enumerate(steps):
            print(step)
            dirname = str(step+1)+"/" # base 1
            t = step*0.1 # fs


            for i_k_point in range(n_k_points):

                chem_pot, evs, evecs, ham1, ovlp1 = emp.parse_timestep_data_kpoint(dirname,i_k_point,i_atom,i_cart,i_spin,parse_evecs=False)
                min_bounds,max_bounds = emp.find_evs_sensible_bounds(evs,chem_pot,max_energy_from_fermi)
                #n_states = np.shape(evecs)[1]


                nac = ham1 - (ovlp1.multiply(chem_pot))
                nac = nac.multiply(1/Bohr)

                for j,back_step in enumerate(steps):

                    if back_step > step:
                        break

                    back_dirname = str(back_step+1)+"/" # base 1
                    t_prime = back_step*0.1 # fs
                    tau = t-t_prime #time difference

                    if (tau)>mem_cutoff: #mem integral cutoff
                        break


                    for j_atom in friction_atoms:
                        for j_cart in range(n_cart):
                    
                            # start = time.time() 
                            chem_pot_b, evs_b, evecs_b, ham1_b, ovlp1_b = emp.parse_timestep_data_kpoint(back_dirname,i_k_point,j_atom,j_cart,i_spin,parse_evecs=False)
                            # end = time.time()
                            # print('Time for 1 parse/ s: '+str(end - start))

                            # start = time.time() 
                            min_bounds_b,max_bounds_b = emp.find_evs_sensible_bounds(evs_b,chem_pot_b,max_energy_from_fermi)
                            # end = time.time()
                            # print('Time for 1 boundary find / s: '+str(end - start))





                            # start = time.time() 
                            nac_back = ham1_b - (ovlp1_b.multiply(chem_pot_b))
                            nac_back = nac_back.multiply(1/Bohr)
                            # end = time.time()
                            # print('Time for 1 nac evaluation / s: '+str(end - start))


                            # if (min_bounds[i_k_point]!=min_bounds_b[i_k_point]):
                            #     print('min bounds mismatch')
                            # if (max_bounds[i_k_point]!=max_bounds_b[i_k_point]):
                            #     print('max bounds mismatch')

                            # start = time.time() 
                            for i_state in range(min_bounds[i_k_point],max_bounds[i_k_point]):
                                for j_state in range(min_bounds_b[i_k_point],max_bounds_b[i_k_point]):


                                    e_diff = (evs[i_k_point,i_state]-evs[i_k_point,j_state])

                                    if (e_diff==0):
                                        continue

                                    # Equation 22 , Head-Gordan Tully 1995
                                    friction_kernel[i_atom,i_cart,j_atom,j_cart,i,j] +=  2.*np.real(nac[i_state,j_state]*np.conjugate(nac_back[j_state,i_state])*\
                                        (1/e_diff) * np.cos(one_over_hbar*(e_diff)*(tau)) * k_weights[i_k_point])

                # end = time.time()
                # print('Time for double loop evaluation / s: '+str(end - start))


        for j_atom in friction_atoms:
            for j_cart in range(n_cart):
                np.savetxt('time_steps_i_atom_{:04d}_i_cart_{:01d}_j_atom_{:04d}_j_cart_{:01d}.txt'.format(i_atom,i_cart,j_atom,j_cart),time_steps)
                np.savetxt('friction_kernel_i_atom_{:04d}_i_cart_{:01d}_j_atom_{:04d}_j_cart_{:01d}.txt'.format(i_atom,i_cart,j_atom,j_cart),friction_kernel[i_atom,i_cart,j_atom,j_cart,:,:])






