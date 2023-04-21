from re import L
import numpy as np
from ase.io import read,write






trajectory = read('../0001_no_bottom.traj@:')

friction_atoms = [0,1]
n_cart = 1



for i_atom in friction_atoms:
    for i_cart in range(n_cart):

        time_steps = np.loadtxt('time_steps_atom_{:04d}_cart_{:01d}.txt'.format(i_atom,i_cart))
        time_step = time_steps[1]-time_step[0]
        friction_kernel = np.loadtxt('friction_kernel_atom_{:04d}_cart_{:01d}.txt'.format(i_atom,i_cart))

        friction_force = np.zeros((time_steps))
        friction_work = np.zeros((time_steps))
        friction_velocities = np.zeros(time_steps)

        friction_mass = trajectory[0].get_masses()[i_atom]

        for i, atoms in enumerate(trajectory):

            friction_velocity = atoms.get_velocities()[i_atom,i_cart]
            friction_velocities[i] = friction_velocity

        for t in range(len(time_steps)):
            friction_force[t]=np.trapz(y=friction_kernel[t,:]*friction_velocities[:],x=time_steps)
        #friction_force *= time_step



        friction_work[:] = friction_force[:] * friction_velocities[:]


        np.savetxt(friction_work,'friction_work_atom_{:04d}_cart_{:01d}.txt'.format(i_atom,i_cart))






