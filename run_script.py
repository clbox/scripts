from ase import Atoms,Atom
from ase.io import read,write
import os
import numpy as np
from ase.io.trajectory import TrajectoryWriter
import sys
from ase.io.trajectory import TrajectoryReader
from ase.db import connect
from ase.calculators.aims import Aims
import os
import time




time_start = time.perf_counter()

max_time = 24*60*60 #24 hours in seconds

time_buffer = 2*60*60 #leave 2 hour buffer
#Want to:
#Run many instances of this script on compute nodes
cwd = os.getcwd()

base_local_dir = cwd + '/calcs/'
species_dir = cwd+'/custom_basis/'
friction_atom_indices  = [0,1]
run_command = 'srun /work/E635/E635/msrvhs2/bin/aims.201122.scalapack.mpi.x > aims.out'
calc1 = Aims(
        aims_command=run_command,
        xc='pbe',
        #output=['k_point_list','friction_matrices'],#,'friction_eigenvectors'],
        sc_accuracy_etot = 1e-6,
        sc_accuracy_rho = 1e-5,
        sc_accuracy_eev = 1e-3,
        #sc_accuracy_forces = 0.001,
        sc_iter_limit = 250,
        charge_mix_param = 0.2,
        occupation_type = ['gaussian', 0.1],
        relativistic = ['atomic_zora','scalar'],
        spin = 'none',
        calculate_all_eigenstates  = True,
        calculate_friction = 'numerical_friction',
        friction_numeric_disp = 0.0025,
        friction_broadening_width = 0.6,
        friction_iter_limit = 90,
        friction_temperature = 300,
        friction_output_couplings = True,
        KS_method = 'parallel',
        k_grid = [9,9,1],
        species_dir = species_dir,
        friction_atoms = friction_atom_indices,
       # restart_aims ='wvfn.dat'
)


class cd: 
    """Context manager for changing the current working directory"""
    """https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python/13197763#13197763"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def write_calc(row,calc):
    import time
    #attaches specified calculator to specified atoms object. Creates control.in, geometry.in and submit script
    #for specified cluster and sends to that cluster and submits the job
    idir = str(row.id)+'/'

    atoms = con.get_atoms(id=row.id)

    #make local subdirectory if doesnt exist
    if not os.path.exists(base_local_dir+idir):
        os.mkdir(base_local_dir+idir)


    with cd(base_local_dir+idir):
        try:
            calc.write_input(atoms)
            return
        except:
            print('error writing, setting calc status to -1')
            con.update(row.id, calculation_status=-1)
            return

def calc_COM(atoms,friction_atom_indices):
    friction_masses = atoms.get_masses()[friction_atom_indices]

    friction_positions = atoms.get_positions()[friction_atom_indices]


    COM = ((friction_positions[0,:] * friction_masses[0]) + (friction_positions[1,:] * friction_masses[1])) / \
        (friction_masses[0]+friction_masses[1])


    return COM

def run_calc(row,calc):
    idir = str(row.id)+'/'

    atoms = con.get_atoms(id=row.id)
    #make local subdirectory if doesnt exist
    if not os.path.exists(base_local_dir+idir):
        os.mkdir(base_local_dir+idir)

    with cd(base_local_dir+idir):
        friction_tensor = calc.get_friction_tensor(atoms)

    return friction_tensor






con = connect('traj_database.db')

#with connect('traj_database.db') as con:
for i in range(1,len(con)+1):


    #Check if enough time is left to run the calculation
    print('Checking elapsed time')
    time_now = time.perf_counter()
    time_elapsed = time_now - time_start
    if max_time - time_elapsed < time_buffer:
        print('Exiting, not enough walltime left')
        break




    print('Checking entry ', i)
    row = con.get(id=i)

    if row.calculation_status == 0:
        con.update(id=i,calculation_status=1)
        #Dont want to be connected to database whole time running calc
        #Ideally would exit this loop
        print('Preparing calculation')

        #atoms = con.get_atoms(id=row.id)
        #COM = calc_COM(atoms,friction_atom_indices)
        # if COM[2] > 4.0:
        #     print('COM > 4.0')
        #     calc = calc2
        #     magmoms = np.zeros((len(atoms)))
        #     magmoms[1] = 1.0
        #     atoms.set_initial_magnetic_moments(magmoms)
        #     #con.update(id=i,magmoms=magmoms)
        # else:
        #     print('COM < 4.0')
        #     calc = calc1
        calc = calc1
        #con.update(id=i,atoms=atoms)
        #row = con.get(id=i)
        #write_calc(row,calc)

        try:
            friction_tensor = run_calc(row,calc)
        except:
            print('error writing, setting calc status to -1')
            con.update(id=row.id, calculation_status=-1)
            continue


        con.update(id=i,calculation_status=2)
        continue


    elif row.calculation_status == 1:
        print('Checking calculation progress')

        #Query progress


    elif row.calculation_status == 2:
        print('This calculation has finished')
        continue


    elif row.calculation_status == -1:
        print('This calculation has an error')
        continue  



        


        