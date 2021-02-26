import gamma.friction_coupling as fc
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

aims1 = 'aims.out'
gamma_files1 = glob.glob("friction_gamma*")


a = fc.friction_gamma_parser(aims_file=aims1,gamma_files=gamma_files1)

#########CALC friction_tensor.out#################
b = fc.friction_tensor(a.ks,a.coords,a.eis,a.ejs,a.couplings,a.kweights,a.chem_pot,a.friction_masses,300,0.6,nspin=1)
tensor = b.calc_tensor()
print(tensor)

