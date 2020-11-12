import gamma.friction_coupling as fc
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

aims1 = 'reactant/aims.out'
aims2 = 'adsorption/aims.out'
gamma_files1 = glob.glob("reactant/gamma-couplings*")
gamma_files2 = glob.glob("adsorption/gamma-couplings*")


b = fc.friction_gamma_parser(aims_file=aims1,gamma_files=gamma_files1)
a = fc.friction_gamma_parser(aims_file=aims2,gamma_files=gamma_files2)



# c = fc.calc_time_tensor(a.ks,a.coords,a.eis,a.ejs,a.couplings,a.kweights,a.chem_pot,friction_masses=a.friction_masses,
#     temp=300,sigma=0.01,nspin=1,min_e=-8,max_e=-2,npoints=100)

# tensor = c.evaluate_tensor()
# print(tensor)

d = fc.calc_time_tensor(a.ks,a.coords,a.eis,a.ejs,a.couplings,a.kweights,a.chem_pot,friction_masses=a.friction_masses,
    temp=300,sigma=0.01,nspin=1,min_e=0,max_e=2.4,npoints=100)
tensor = d.evaluate_tensor2()
print(tensor)

x = np.linspace(0,1,100)
y = fc.gaussian_function(x,0,0.01)
y2 = y/(np.sum(y)*(x[1]-x[0]))
fig,ax = plt.subplots(1,1)
ax.plot(x,y)
ax.plot(x,y2)
ax.plot(x,fc.gaussian_norm(x,0.01))
fig.savefig('gauss.pdf')
# print('Max excitation')
# print(np.max(a.ejs-a.eis))
# print('Max imaginary component of A')
# print(np.max(np.imag(A)))
# fig,ax = plt.subplots(1,1)
# pos = ax.imshow(np.real(A[5,5,:,:]),cmap='jet',interpolation='none',extent=(-8,-2,-2,-8))
# fig.colorbar(pos, ax=ax)
# ax.set_ylabel(r'$\epsilon_e$ / eV')
# ax.set_xlabel(r'$\epsilon_h$ / eV')
# fig.savefig('heat.pdf')
#########CALC friction_tensor.out#################
#b = fc.friction_tensor(a.ks,a.coords,a.eis,a.ejs,a.couplings,a.kweights,a.chem_pot,a.friction_masses,300,0.6,nspin=1)
#tensor = b.calc_tensor()
#print(tensor)

#############Discretize to Gamma####################
#c = fc.calc_gamma(a.ks,a.coords,a.eis,a.ejs,a.couplings,a.kweights,a.chem_pot,sigma=0.01,nspin=1,min_e=-8,max_e=-1,npoints=10)
#gamma = c.dos_binning()


###########Calc Lambda####################
