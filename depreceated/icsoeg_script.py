
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from scipy.interpolate import interp1d,UnivariateSpline
import scipy.integrate as integrate
import scipy.linalg as linalg
import pandas as pd


# In[2]:


#https://codereview.stackexchange.com/questions/82010/averaging-lists-of-values-with-duplicate-keys
def avg_dups(eigenvalues, values):
    folded, indices, counts = np.unique(eigenvalues, return_inverse=True, return_counts=True)
    output = np.zeros((folded.shape[0]))
    np.add.at(output, indices, values)
    output /= counts[:]

    return folded, output

boltzmann_kB = 3.166811429e-6 # hartree per kelvin (wikipedia)
def fermi_derivative(x,x0,T):
        fermi_derivative = np.exp((x-x0)/(boltzmann_kB*T))
        fermi_derivative = fermi_derivative / (((fermi_derivative+1)**2)*T*boltzmann_kB)
        
        return fermi_derivative
hbar = 1 
hartree = 27.2113845 #* to convert hartree to eV
ps = 2.418884326505e-5 #* to convert ps to a.u
max_energy = 5.4/hartree
x = np.linspace(-max_energy,max_energy,1000)
friction_atoms = np.array([16,17])

# In[3]:



filenames = []
kpoint = 0
colours = ("#009E73","#0072B2", "#D55E00", "#CC79A7","#E69F00","#F0E442","#000000","#56B4E9")
cwd = os.getcwd()

filenames = []
for filename in os.listdir(cwd):
    if "nacs-spectrum" in filename:
        if filename.endswith(".out"):
            filenames.append(filename)
        else:
            continue
    else:
        continue

filenames.sort()
with open(filenames[0], "r") as f:
    for line in f:
        if "Friction" in line:
            dimension = int(line.split()[3])

print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
elements = (((dimension*dimension)-dimension)/2)+dimension
print("There are " + str(elements) + " coupling components")

component_dict = dict()
header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point"]
for filename in filenames[1:]: #SKIPPING first k point
    kpoint = int(filename[15:-4])
    with open(cwd+'/'+filename, "r") as f:
        file_data = []
        for line in f:
            if any(x in line for x in header):
                continue
            elif "Friction" in line:
                component = int(line.split()[3]) 
            else:
                if component in component_dict:
                    component_dict[component].append(line)
                else:
                    component_dict[component] = [line]


# In[ ]:


#Method Dos: Avg dups, cumulative sum, fit, couple then calculate the friction
#2.1 
#CALCULATE TENSOR AND VISUALIZE FIT
dfde = fermi_derivative(x,0,300)
tensor = np.zeros((6,6))
fits_r = []
fits_i = []
for i in range(dimension):  
    a = np.loadtxt(component_dict[i+1])
    e,real = avg_dups(a[:,0],a[:,1])
    e,imag = avg_dups(a[:,0],a[:,2])
    e_inds = e.argsort()
    e_s= (e)[e_inds]
    A = (real)[e_inds]+1j*(imag)[e_inds] 
    #cumsum fit & derivative

    sum_r = np.poly1d(np.polyfit(e_s, np.cumsum(A),4))
    sum_fr = np.polyder(sum_r)     
    
    fits_r.append(sum_fr)
    #fits_i.append(sum_fi)
            
for i in range(dimension):
    for j in range(dimension):
        if i <= j:               
            ##CUMULATIVE SUM
            #y = (fits_r[i](x)-1j*fits_i[i](x))*(fits_r[j](x)+1j*fits_i[j](x))*dfde
            y = np.conj(fits_r[i](x))*(fits_r[j](x))*dfde
            tensor[i,j] = integrate.trapz(np.real(y),x)*hbar*np.pi
            tensor[j,i] = np.conj(tensor[i,j])


# In[ ]:


hbar_to_js = 1.05457180013e-34
a0_to_m = 5.29177e-11
au = 1.6605e-27 #kg

tensor_si = tensor * hbar_to_js #hbar to Js
tensor_si /= (a0_to_m)**2

row_labels = ['O_x', 'O_y', 'O_z', 'C_x','C_y','C_z']
column_labels = ['O_x', 'O_y', 'O_z','C_x','C_y','C_z']

coord = [16,16,16,12,12,12]

for i in range(dimension):
    for j in range(dimension):
            tensor_si[i,j] /= np.sqrt(coord[i]*au*coord[j]*au)
            
tensor_si *= 1e-12 #s^-1 to ps^-1
df_si = pd.DataFrame(tensor_si,columns=column_labels, index=row_labels)
df_si.round(2)

eigvals,eigvecs = np.linalg.eig(tensor_si)
print('Friction mode lifetimes / ps')
for i,ii in enumerate(eigvals):
    print('{:1f}').format(np.real(1/eigvals[i]))
    
with open('friction_tensor_icsoeg.out','w') as f:
    f.write(' # friction_tensor in 1/ps \n')
    i=0 
    for n in friction_atoms:
        for i_cart in range(3):
            i += 1
            f.write("'# n_atom  {0:2d} i_cart  {1:2d}  index  {2:2d}\n".format(n,i_cart+1,i))
            for j in range(dimension):
                f.write(" {:12.6f} ".format(np.real(tensor[i-1,j])))
            f.write('\n')


# In[ ]:


#Method tres: Avgs dups, Cumulative sum, fit then coupling per k point. Average friction over all k points.

