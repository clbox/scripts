from ase.build import molecule, bulk, surface, add_adsorbate, fcc111, hcp0001
from ase.io import read,write
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ase.calculators.aims import Aims
import numpy as np
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState

matplotlib.use('PDF')
matplotlib.rcParams['font.sans-serif'] = "Arial"

# User parameters
element = 'Pt'
structure = 'fcc'
approx_constant = 3.6 # Check Haas 2009 paper for example
number_of_values = 5
aims_dir = '/home/chem/msrvhs/software/aims/FHIaims/'
species_dir = aims_dir+'species_defaults/defaults_2020/'+'tight/'
run_command = 'srun '+aims_dir+'build_master/aims.210825.scalapack.mpi.x > aims.out'


# Set up calculator
calc1 = Aims(
        aims_command=run_command,
        xc='rpbe',
        sc_accuracy_etot = 1e-6,
        sc_accuracy_rho = 1e-5,
        sc_accuracy_eev = 1e-3,
        #sc_accuracy_forces = 0.001,
        sc_iter_limit = 250,
        charge_mix_param = 0.05,
        occupation_type = ['gaussian', 0.1],
        relativistic = ['atomic_zora','scalar'],
        spin = 'none',
        KS_method = 'parallel',
        k_grid = [3,3,3],
        species_dir = species_dir,
)

# Calculate raw data
atoms = bulk(element, structure, a=approx_constant)
atoms.calc = calc1
cell = atoms.get_cell()

traj = Trajectory('BM.traj', 'w')
for x in np.linspace(0.95, 1.05, number_of_values):
    atoms.set_cell(cell * x, scale_atoms=True)
    atoms.get_potential_energy()
    traj.write(atoms)



# Extract volumes and energies:
configs = read('BM.traj')  # read 5 configurations
volumes = [atoms.get_volume() for atoms in configs]
energies = [atoms.get_potential_energy() for atoms in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')

lattice_constant = (4 * v0)**(1 / 3.0)
print('Optimised lattice constant is '+str(lattice_constant)+' / Angstrom')


eos.plot('eos.png')


