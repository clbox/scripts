import gamma.friction_coupling as fc
from scipy.integrate import simps
import glob
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.style.use('clb_publication')
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


dir_name1 = '300/'
dir_name2 = '301/'


gamma_files1 = glob.glob(dir_name1+'/*friction_gamma*.out')
gamma_files2 = glob.glob(dir_name2+'/*friction_gamma*.out')

parser = fc.friction_output_parser_2021()
eigenvalues1, couplings1, n_excits1, k_points1 = parser.parse_gamma_couplings(gamma_files1)
eigenvalues2, couplings2, n_excits2, k_points2 = parser.parse_gamma_couplings(gamma_files2)

chem_pot1 = parser.parse_chem_pot(dir_name1+'aims.out')
chem_pot2 = parser.parse_chem_pot(dir_name2+'aims.out')
print('Chemical potentials / eV : '+ str(chem_pot1) + ' ' + str(chem_pot2) + ' Difference / eV : ' +str(chem_pot1-chem_pot2))


fig, ax = plt.subplots(1, 1)

# print(np.shape(eigenvalues2[0]))
# print(eigenvalues2[0][:,0])
for eig in range(len(eigenvalues1[0][:,0])):
    print(np.array((eigenvalues1[0][eig,0],eigenvalues2[0][eig,0])))
    ax.plot(np.array((eigenvalues1[0][eig,0],eigenvalues2[0][eig,0])),color='dodgerblue', marker=" ",linestyle='-')

for eig in range(len(eigenvalues1[0][:,1])):
    ax.plot(np.array((eigenvalues1[0][eig,1],eigenvalues2[0][eig,1])),color='orangered', marker=" ",linestyle='-')   

ax.plot(np.array((chem_pot1,chem_pot2)), color='black', marker=" ",linestyle='-')

# ax.hist(eigenvalues1[0][:,0],bins=bins,color='orangered')
# ax.hist(eigenvalues2[0][:,0],bins=bins,color='dodgerblue')
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(10))

ax.set_ylabel(r'$\epsilon$ / $\mathrm{eV}$')
ax.set_xlabel(r'$\mathrm{Structure}$')
# ax.set_ylabel(r'$\Lambda_{0,0}(\epsilon)$ / ')
fig.savefig('eigenvalues_300_301.pdf',transparent=False,bbox_inches='tight')


# ...... Energy tensor
e_cutoff = 1.0 #eV
n_points = 200
time_cutoff = 10 #fs
n_time_points = 200
temp = 300 # K
sigma = 0.1 # eV

#time_cutoffs = [0.1,3,9,12,15,50,100]

fig, ax = plt.subplots(1, 1)
sigmas = [0.01,0.05,0.1,0.3]
for sigma in sigmas:
    friction_masses = parser.parse_fiction_masses(dir_name1+'aims.out')
    ctt = fc.calc_time_tensor_2021(chem_pot1,e_cutoff,friction_masses,temp=temp,sigma=sigma)
    ex_energy_grid, ex_energy_grid_tensor = ctt.calculate_ex_energy_grid_tensor(eigenvalues1,eigenvalues2,couplings1,couplings2,n_excits1,n_excits2,k_points1[:,3],n_points)
    ax.plot(ex_energy_grid,ex_energy_grid_tensor[0,0,:],marker='None',label=str(sigma))

# ...... Time tensor

#time_axis, time_tensor = ctt.evaluate_time_tensor(ex_energy_grid,ex_energy_grid_tensor,time_cutoff,n_time_points)
# plot

#ax.plot(time_axis,time_tensor[0,0,:],label=str(e_cutoff))#,color='orangered')

        
ax.legend()
ax.set_xlim(0,0.8)
ax.set_ylim(-0.2,0.2)

ax.set_xlabel(r'$\epsilon$ / $\mathrm{eV}$')
ax.set_ylabel(r'$\Lambda_{0,0}(\epsilon)$ / ')

fig.savefig('energy_tensor_300_301.pdf',transparent=False,bbox_inches='tight')