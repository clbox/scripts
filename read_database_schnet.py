from ase.io import read,write
from ase.visualize import view
import numpy as np
from ase import Atoms
import os
from ase.db import connect
from ase.calculators.aims import Aims
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from numpy import format_float_scientific
import pandas as pd
import schnetpack as spk

def read_raw_spectra(path):
    filenames = []
    header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point"]
    for filename in os.listdir(path):
        if "nacs-spectrum" in filename:
            if filename.endswith(".out"):
                filenames.append(filename)
            else:
                continue
        else:
            continue

    filenames.sort()

    #open a file and get dimension
    head_count=0
    with open(path+filenames[0], "r") as f:
        for line in f:
            if "Friction" in line:
                dimension = int(line.split()[3])
                head_count += 1
            if any(x in line for x in header):
                continue
    #get spin status
    elements = (((dimension*dimension)-dimension)/2)+dimension
    if elements < head_count:
        n_spin = 2


    component_dict = dict()
    for filename in filenames:
        kpoint = int(filename[15:-4])
        with open(path+filename, "r") as f:
            file_data = []
            for line in f:
                if any(x in line for x in header):
                    continue
                elif "Friction" in line:
                    component = line.split()[3]+line.split()[4]
                else:
                    if component in component_dict:
                        component_dict[component].append(line)
                    else:
                        component_dict[component] = [line]
    
    return component_dict,dimension

def read_memory_kernel(path):

    head_count =0
    header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction"] #skip lines
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                dimension = int(line.split()[3])
                head_count += 1
            if "Discretization" in line:
                discretization=float(line.split()[-1])
            if any(x in line for x in header):
                continue
            max_e = float(line.split()[0])
    elements = (((dimension*dimension)-dimension)/2)+dimension
    if elements < head_count:
        n_spin = 2

    bins=np.zeros((int(max_e/discretization)+1))
    re_memory_kernel=np.zeros((dimension,dimension,len(bins)))
    im_memory_kernel=np.zeros_like(re_memory_kernel)
    
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                i = int(line.split()[3])
                j = int(line.split()[4])
                head_count += 1
                c=0
            if any(x in line for x in header):
                continue
            else:
                re_memory_kernel[i-1,j-1,c]=float(line.split()[1])
                im_memory_kernel[i-1,j-1,c]=float(line.split()[2])
                bins[c]=float(line.split()[0])
                c +=1
    return(bins,re_memory_kernel,im_memory_kernel,dimension,max_e)

def string2array(string):
    """
    
    """
    dimension = len(string.split(']\n'))
    onedarray = np.fromstring((string.replace('[',' ').replace(']\n',' ')),dtype=float,sep=' ')
    return onedarray.reshape(dimension,dimension)


# PARSING DATA 
n_geoms = 797
available_properties = ["raw_coupling_energies","raw_coupling_elements","smear_frequency_energies","smeared_frequency_friction","markov_friction_tensor",]

dataset = spk.data.AtomsData("CO_Cu100_FT_ML.db", available_properties=available_properties)
raw_dataset = connect('raw_database.db')

for idx in range(1,n_geoms+1):
#for idx in range(1,2):
    print(idx)
    property_list = {}
    mol = raw_dataset.get_atoms(id=idx)
    
    #Markov
    row = raw_dataset.get(id=idx)
    try:
        markov_tensor = string2array(row.get('ft1'))
        property_list.update({"markov_friction_tensor": markov_tensor})
    except:
        print('no markov')
        continue
    #RAW

    try:
        raw_coupling_data,dimension = read_raw_spectra('ml/'+str(idx)+'/')
        raw_coupling_energies = np.zeros((len(raw_coupling_data['11'])))
        raw_tensor = np.zeros((dimension,dimension,len(raw_coupling_data['11'])))
        for i in range(dimension):
            for j in range(dimension):
                if j>=i:
                    a = np.loadtxt(raw_coupling_data[str(i+1)+str(j+1)])
                    raw_coupling_energies[:],raw_tensor[i,j,:] = a[:,0],a[:,1]
                    raw_tensor[j,i,:]=raw_tensor[i,j,:]
        property_list.update({"raw_coupling_energies": raw_coupling_energies,"raw_coupling_elements": raw_tensor})
    except:
        print('no raw')
        continue
    #Smeared
    try:
        bins,re_memory_kernel,im_memory_kernel,dimension,max_e = read_memory_kernel('ml/'+str(idx)+'/friction_memory_kernel.out') 
        for i in range(dimension):
            for j in range(dimension):
                if j>=i:
                    re_memory_kernel[j,i]=re_memory_kernel[i,j]
        property_list.update({"smear_frequency_energies": bins,"smeared_frequency_friction": re_memory_kernel})
    except:
        print('no smear')
        continue

    dataset.add_system(mol, **property_list)
    print('added')
