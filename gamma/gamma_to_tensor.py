import gamma.friction_coupling as fc
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

aims1 = 'aims.out'
gamma_files1 = glob.glob("friction_couplings*")


a = fc.friction_gamma_parser(aims_file=aims1,gamma_files=gamma_files1)

print(len(a.couplings))
print(len(a.coords))
print(len(a.eis))
print(len(a.ejs))
print(len(a.ks))

print(np.where((a.coords == 0)))
print(np.where((a.coords == 1)))
i_idx = np.where((a.coords == 0) & (a.ks == 0))[0]
j_idx = np.where((a.coords == 1) & (a.ks == 0))[0]

print('---------')
print(len(a.eis[i_idx]))
print(len(a.ejs[j_idx]))
print('---------')
#########CALC friction_tensor.out#################
b = fc.friction_tensor(a,300,0.6,nspin=1)
tensor = b.calc_tensor()
print(tensor)

