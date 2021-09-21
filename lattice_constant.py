from ase.build import molecule, bulk, surface, add_adsorbate, fcc111, hcp0001
from ase.io import read,write
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

matplotlib.use('PDF')
matplotlib.rcParams['font.sans-serif'] = "Arial"

# User parameters
element = 'Ru'
structure = 'hcp'
min_constant = 2.5
max_constant = 2.7
number_of_values = 30
aims_dir = '/home/chem/msrvhs/software/aims/FHIaims/'
species_dir = aims_dir+'species_defaults/'+'tight/'
run_command = 'srun '+aims_dir+'build_dfpt2/aims.210513.scalapack.mpi.x > aims.out'

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
constants_to_test = np.linspace(min_constant,max_constant,number_of_values)
energies = np.zeros_like(constants_to_test)
for i, constant in enumerate(constants_to_test):
    atoms = bulk(element, structure, a=constant)
    atoms.calc = calc1

    energy = atoms.get_potential_energy()

    print(str(constant)+' '+str(energy))

    energies[i] = energy


# Plot raw data
fig, ax = plt.subplots(1, 1)
ax.plot(constants_to_test,energies,linestyle='none',marker='o',mfc='none',color='black',label='SCF energy')

# Carry out fit and plot
def x_squared(x,a,c):
    y = a * x^2 + c
    return y 
popt, pcov = curve_fit(x_squared, constants_to_test, energies)
fit_x = np.linespace(2,5,2.7,300)
fit_data = x_squared(fit_x,popt)
ax.plot(fit_x,fit_data,linestyle='--',color='black',label='Fit')

# Find minimum energy lattice constant
idx = np.argmin(fit_data)
lattice_constant = fit_x[idx]
print('Optimised lattice constant is '+str(lattice_constant)+' / Angstrom')

# Output plot
ax.set_xlabel(r'Lattice constant / $\mathrm{\AA{}}$')
ax.set_ylabel('Potential energy / eV')

ax.legend(fancybox=True,framealpha=0)

fig.set_figheight(3.5)
fig.set_figwidth(4.)
fig.savefig('lattice_constant.pdf',transparent=False,bbox_inches='tight')