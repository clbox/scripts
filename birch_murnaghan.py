from ase.build import molecule, bulk, surface, add_adsorbate, fcc111, hcp0001
from ase.io import read,write
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ase.eos import calculate_eos

matplotlib.use('PDF')
matplotlib.rcParams['font.sans-serif'] = "Arial"

# User parameters
element = 'Ru'
structure = 'hcp'
start_constant = 2.8
aims_dir = '/home/chem/msrvhs/software/aims/FHIaims/'
species_dir = aims_dir+'species_defaults/'+'tight/'
run_command = 'srun '+aims_dir+'build_dfpt2/aims.210513.scalapack.mpi.x > aims.out'

bulk_atoms = bulk(element, structure, a=start_constant)
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

bulk_atoms.calc = calc1
eos = calculate_eos(bulk_atoms, trajectory='bm.traj')
v, e, B = eos.fit()
a = (4 * v)**(1 / 3.0)
print('{0:.6f}'.format(a))

