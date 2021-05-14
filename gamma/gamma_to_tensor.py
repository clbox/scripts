import gamma.friction_coupling as fc
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

aims1 = 'aims.out'
gamma_files1 = glob.glob("*gamma*k*out")


a = fc.friction_gamma_parser(aims_file=aims1,gamma_files=gamma_files1)


#########CALC friction_tensor.out#################
b = fc.friction_tensor(a,300,0.6,nspin=1)
tensor = b.calc_tensor()
print(tensor)

modes = ['default','double_delta','double_delta_half_sigma','prb_print','no_norm']
for mode in modes:
    print('---------'+mode+'-----------')
    tensor = b.calc_tensor()
    print(tensor)
    fl = open('friction_tensor_'+mode+'.out', 'w')
    for i in range(np.shape(tensor)[0]):
        string = ''
        for j in range(np.shape(tensor)[1]):
            string = '{0:.15f}'.format(tensor[i,j])
            fl.write(string+'\n')
