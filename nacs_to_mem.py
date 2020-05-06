#!/usr/bin/env python
from ase.io import read,write
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from numpy import format_float_scientific
import pandas as pd
import sys

def read_spectral_data(path):


    #maximum energy of data in nacs_spectrum_kxxx.out
    #Grid density to discretize the fit on
    #Number of knots of the piecewise cubic spline fit
    #eigenvalue mode investigating starting from 0. E.g for 6x6 tensor it will be 0 or 1 or 2 .. 5 
    #Choose the broadening or range of broadenings to investigate
    # # Parser and calculation
    #frequenccy evaluated at in eV
    n_spin = 1
    head_count = 0

    header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction"] #skip lines                    
    ##########################################################################
    filenames = []
    for filename in os.listdir(path):
        if "nacs-spectrum" in filename:
            if filename.endswith(".out"):
                filenames.append(filename)
            else:
                continue
        else:
            continue

    filenames.sort()
        
    with open(path+filenames[0], "r") as f:
        for line in f:
            if "Friction" in line:
                dimension = int(line.split()[3])
                head_count += 1
            if any(x in line for x in header):
                continue
            max_e = float(line.split()[0])
    #print("Friction max energy = "+str(max_e))
    #print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
    elements = (((dimension*dimension)-dimension)/2)+dimension
   # print("There are " + str(elements) + " coupling components")
    if elements < head_count:
        n_spin = 2
        #print("This system is spin unrestricted")
    zpf = np.zeros((dimension,dimension),dtype=complex)
    file_counter = -1

    for filename in filenames: #FILENAME LOOP#######
        file_counter += 1

        with open(path+filename, "r") as f:
            file_data = []
            for line in f:
                if any(x in line for x in header):
                    continue
                else:
                    file_data.append(line)

        a = np.loadtxt(file_data)

        if file_counter == 0:
            data_sum = np.zeros([len(a[:,0])//n_spin,3])
        data_split = np.vsplit(a,elements*n_spin)
        c = 0
        i = 0
        for chunk in data_split:
            c += 1
            if c > elements:
                i = 0
                c = 0
            for row in chunk:
                data_sum[i,0] = row[0]
                data_sum[i,1] += row[1]
                data_sum[i,2] += row[2]
                i += 1

    ##END FILENAME LOOP ######################
    data_split = np.vsplit(data_sum,elements)   
    return(data_split,dimension,max_e)
         


paths = sys.argv[1:]
for path in paths:
    #try:
    data_split,dimension,max_e = read_spectral_data(path)
    #except:
    #    print('Couldnt  read data - continuing')
    #    continue
    print('max energy = {} eV'.format(max_e))
    c = -1
    i_cart = -1
    i_atom = 0
    for i in range(dimension):
        i_cart = i_cart + 1
        if (i_cart>2):
            i_cart = 0
            i_atom = i_atom + 1
        j_atom = 0
        j_cart = -1          
        for j in range(dimension):
            j_cart = j_cart + 1
            if (j_cart>2):
                j_cart = 0
                j_atom = j_atom + 1
                
            if j >= i:
                c +=1         
                n_bins = len((data_split[c])[:,0])
                e_axis = (data_split[c])[:,0]
                real = (data_split[c])[:,1]
                imag = (data_split[c])[:,2]
    
                with open(path+'friction_memory_kernel.out','a+') as f:
                    if i==0 and j==0:
                        f.write("No of components        {}\nDiscretization length in eV   {}\nNumber of Bins      {}\n"\
                                "Excitation energy in eV   Re(Coupling element)   Im(Coupling element) in 1/ps\n" \
                                "============================================================================\n".format(dimension,max_e/(n_bins-1),n_bins))
                            
                    f.write('Friction component for '+str(i+1)+' '+str(j+1)+'\n')
                    for e_counter in range(n_bins):
                        f.write('{0:1.6E}    {1:1.6E}  {2:1.6E}\n'.format(e_axis[e_counter],real[e_counter],imag[e_counter]))
