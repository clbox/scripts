import gamma.friction_coupling as fc
import glob
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.style.use('clb_publication')
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


dir_name1 = '1/'
dir_name2 = '2/'


gamma_files1 = glob.glob(dir_name1+'/*friction_gamma*.out')
gamma_files2 = glob.glob(dir_name2+'/*friction_gamma*.out')

parser = fc.friction_output_parser_2021()
eigenvalues1, couplings1, n_excits1, k_points1 = parser.parse_gamma_couplings(gamma_files1)
eigenvalues2, couplings2, n_excits2, k_points2 = parser.parse_gamma_couplings(gamma_files2)

chem_pot1 = parser.parse_chem_pot(dir_name1+'aims.out')
chem_pot2 = parser.parse_chem_pot(dir_name2+'aims.out')
print('Chemical potentials / eV : '+ str(chem_pot1) + ' ' + str(chem_pot2) + ' Difference / eV : ' +str(chem_pot1-chem_pot2))



# ...... Energy tensor
e_cutoff = 0.1 #eV
n_points = 100

friction_masses = parser.parse_fiction_masses(dir_name1+'aims.out')
ctt = fc.calc_time_tensor_2021(chem_pot1,e_cutoff,friction_masses,temp=300,sigma=0.01)
ex_energy_grid, ex_energy_grid_tensor = ctt.calculate_ex_energy_grid_tensor(eigenvalues1,eigenvalues2,couplings1,couplings2,n_excits1,n_excits2,k_points1[:,3],n_points)

fig, ax = plt.subplots(1, 1)

ax.plot(ex_energy_grid,ex_energy_grid_tensor[0,0,:])#,color='orangered')

ax.set_xlim(0,e_cutoff)

ax.set_xlabel(r'$\epsilon$ / $\mathrm{eV}$')
ax.set_ylabel(r'$\Lambda_{0,0}(\epsilon)$ / ')

# fig.set_figheight(2.5)
# fig.set_figwidth(2.5)

fig.savefig('energy_tensor.pdf',transparent=False,bbox_inches='tight')



# ...... Time tensor
time_cutoff = 10 #fs
n_time_points = 100
time_axis, time_tensor = ctt.evaluate_time_tensor(ex_energy_grid,ex_energy_grid_tensor,time_cutoff,n_time_points)

# plot
fig, ax = plt.subplots(1, 1)

ax.plot(time_axis,time_tensor[0,0,:])#,color='orangered')

ax.set_xlim(0,10)

ax.set_xlabel(r'$\tau$ / $\mathrm{fs}$')
ax.set_ylabel(r'$\Lambda_{0,0}(\tau)$ / ')

# fig.set_figheight(2.5)
# fig.set_figwidth(2.5)
fig.savefig('time_tensor.pdf',transparent=False,bbox_inches='tight')