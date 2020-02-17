import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from scipy.interpolate import interp1d,UnivariateSpline,CubicSpline,Rbf
import scipy.integrate as integrate
import scipy.linalg as linalg
import pandas as pd
from basis_expansions import (Binner, Polynomial, 
                              LinearSpline, CubicSpline,
                              NaturalCubicSpline)
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from coolvib.routines.output import print_matrix
import sys
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline


def get_natural_cubic_spline_model(x, y, minval=None, maxval=None, n_knots=None, knots=None):
    """
    Get a natural cubic spline model for the data.

    For the knots, give (a) `knots` (as an array) or (b) minval, maxval and n_knots.

    If the knots are not directly specified, the resulting knots are equally
    space within the *interior* of (max, min).  That is, the endpoints are
    *not* included as knots.

    Parameters
    ----------
    x: np.array of float
        The input data
    y: np.array of float
        The outpur data
    minval: float 
        Minimum of interval containing the knots.
    maxval: float 
        Maximum of the interval containing the knots.
    n_knots: positive integer 
        The number of knots to create.
    knots: array or list of floats 
        The knots.

    Returns
    --------
    model: a model object
        The returned model will have following method:
        - predict(x):
            x is a numpy array. This will return the predicted y-values.
    """

    if knots:
        spline = NaturalCubicSpline(knots=knots)
    else:
        spline = NaturalCubicSpline(max=maxval, min=minval, n_knots=n_knots)

    p = Pipeline([
        ('nat_cubic', spline),
        ('regression', LinearRegression(fit_intercept=True))
    ])

    p.fit(x, y)

    return p

def get_cubic_spline_model(x, y, minval=None, maxval=None, n_knots=None, knots=None):
    """
    Get a cubic spline model for the data.

    For the knots, give (a) `knots` (as an array) or (b) minval, maxval and n_knots.

    If the knots are not directly specified, the resulting knots are equally
    space within the *interior* of (max, min).  That is, the endpoints are
    *not* included as knots.

    Parameters
    ----------
    x: np.array of float
        The input data
    y: np.array of float
        The outpur data
    minval: float 
        Minimum of interval containing the knots.
    maxval: float 
        Maximum of the interval containing the knots.
    n_knots: positive integer 
        The number of knots to create.
    knots: array or list of floats 
        The knots.

    Returns
    --------
    model: a model object
        The returned model will have following method:
        - predict(x):
            x is a numpy array. This will return the predicted y-values.
    """

    if knots:
        spline = CubicSpline(knots=knots)
    else:
        spline = CubicSpline(max=maxval, min=minval, n_knots=n_knots)

    p = Pipeline([
        ('nat_cubic', spline),
        ('regression', LinearRegression(fit_intercept=True))
    ])

    p.fit(x, y)

    return p

from math import exp, sin, sqrt, pi
hartree = 27.2113845 #* to convert hartree to eV
ps = 2.418884326505e-5 #* to convert ps to a.u


# # Calculation parameters

# In[4]:

#Choose colours for each k grid
colours = ("#009E73","#0072B2", "#D55E00", "#CC79A7","#E69F00","#F0E442","#000000","#56B4E9")

natural = False
#bin_size of data in nacs_spectrum_kxxx.out
bin_size = 0.01
#maximum energy of data in nacs_spectrum_kxxx.out
#Grid density to discretize the fit on
npoints = 100000
#Number of knots of the piecewise cubic spline fit
nknots = int(sys.argv[2])
#eigenvalue mode investigating starting from 0. E.g for 6x6 tensor it will be 0 or 1 or 2 .. 5 
#Choose the broadening or range of broadenings to investigate
friction_atoms = np.array([17,18]) #like AIMS, starting from 1, for printing only
# # Parser and calculation
freq = float(sys.argv[1]) #frequenccy evaluated at in eV
# In[5]:
n_spin = 1
head_count = 0

header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction"] #skip lines
cwd = os.getcwd()                    
##########################################################################
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
            head_count += 1
        if any(x in line for x in header):
            continue
        max_e = float(line.split()[0])
print("Friction max energy = "+str(max_e))
grid_discretization = max_e/npoints
print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
elements = (((dimension*dimension)-dimension)/2)+dimension
print("There are " + str(elements) + " coupling components")
if elements < head_count:
    n_spin = 2
    print("This system is spin unrestricted")
zpf = np.zeros((dimension,dimension),dtype=complex)
file_counter = -1

for filename in filenames: #FILENAME LOOP#######
    file_counter += 1

    with open(cwd+'/'+filename, "r") as f:
        file_data = []
        for line in f:
            if any(x in line for x in header):
                continue
            else:
                file_data.append(line)
    a = np.loadtxt(file_data)
    if file_counter == 0:
        data_sum = np.zeros([len(a[:,0])/n_spin,3])
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
x =np.linspace(freq,freq+1.0,2)/hartree
    
################################################
i = 0 
j = 0 
c = -1
for i in range(dimension):
    for j in range(dimension):
        if j >= i:
            #print('Friction component for ',i+1,j+1)
            c += 1
            e_axis = (data_split[c])[:,0]/hartree
            real = (data_split[c])[:,1]*ps
            imag = (data_split[c])[:,2]*ps
            if natural:
                model_6 = get_natural_cubic_spline_model(e_axis, real, minval=min(e_axis), maxval=max(e_axis), n_knots=nknots)
            else:
                model_6 = get_cubic_spline_model(e_axis, real, minval=min(e_axis), maxval=max(e_axis), n_knots=nknots)
                
            #NORMAL       
            zpf[i,j] = model_6.predict(x)[0]
            zpf[j,i] = model_6.predict(x)[0]
                        



zpf_relaxs = np.zeros((dimension),dtype=complex)
zpf_lifetimes = np.zeros((dimension))
    
zpf_relaxs[:],eigvecs=np.linalg.eigh(zpf[:,:])
zpf_lifetimes[:] = ps/zpf_relaxs[:]
    
print('Friction mode lifetimes / ps using Interpolation w/o broadening at ' + str(x[0]*hartree)+' eV')
print(zpf_lifetimes)

print('Friction tensor in 1/ps')
i=0
for n in friction_atoms:
    for i_cart in range(3):
        i += 1
        print('# n_atom ', n,' i_cart ', i_cart+1, ' index ', i)

        print(np.real(zpf[i-1,:]/ps))

with open('friction_tensor_'+str(nknots)+'_'+'{:.0f}'.format(freq*1000)+'meV.out','w') as f:
    f.write(' # friction_tensor in 1/ps \n')
    i=0
    for n in friction_atoms:
        for i_cart in range(3):
            i += 1
            f.write("'# n_atom  {0:2d} i_cart  {1:2d}  index  {2:2d}\n".format(n,i_cart+1,i))
            for j in range(dimension):
                f.write(" {:12.6f} ".format(np.real(zpf[i-1,j]/ps)))
            f.write('\n')


