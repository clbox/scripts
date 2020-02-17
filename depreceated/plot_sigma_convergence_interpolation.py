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


#######
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline


# # Helper spline function

# In[2]:


import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
import sys

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


# # Broadening

# In[3]:


from math import exp, sin, sqrt, pi
hartree = 27.2113845 #* to convert hartree to eV
ps = 2.418884326505e-5 #* to convert ps to a.u

def gaussian(x,x0,s):
    try:
        return 1./sqrt(2*pi*s*s)*exp(-0.5*(((x-x0) / (s))**2))
    except OverflowError:
        return 0.0
# Numerical normalisation
def evaluate_delta_function(x_axis, f, x0, sigma, delta_method):
    """
    evaluates the first moment of a function times a delta function
    """
    
    norm = 0.0
    result = 0.0
    for i, x in enumerate(x_axis):
        delta = delta_function(x,x0,sigma, delta_method)
        norm += delta
        result += f[i]*delta
    result/= norm

    return result
def delta_function(x,x0,s, method):
    """
    function generator that returns a function, which 
    yields true or false if its argument is inside the given 
    window defined by minimum and maximum.
    """
    epsilon = 1E-100
    from math import isnan

    if method is 'gaussian':
        dirac_weight = gaussian(x,x0,s)
    elif method is 'square':
        dirac_weight = square(x,x0,s)
    elif method is 'sine':
        dirac_weight = sine(x,x0,s)
    elif method is 'lorentzian':
        dirac_weight = lorentzian(x,x0,s)
    elif method is 'sine_VSB':
        dirac_weight = sine_VSB(x,x0,s)
    elif method is 'squashed_fermi':
        dirac_weight = squashed_fermi(x,x0,s)
    else:
        raise NotImplementedError('delta method {0} is unknown'.format(method))

    if dirac_weight<epsilon:
        dirac_weight = 0.0
    return dirac_weight


# Analytic normalisatiion
def evaluate_delta_function_2(x_axis, f, x0, sigma, delta_method):
    """
    evaluates the first moment of a function times a delta function
    """
    
    result = 0.0
    for i, x in enumerate(x_axis):
        #delta = delta_function(x,x0,sigma, delta_method)/(0.5*(1.-np.math.erf((-x/sigma)*(1./np.sqrt(2.)))))
        #result += f[i]*delta
        result += f[i]*delta_function(x, x0, sigma, delta_method)/            gaussian_norm(x,sigma)
        #print(delta_function(x, x0, sigma, delta_method)/gaussian_norm(x,sigma))
            #(0.5*(1.-np.math.erf((-x/sigma)*(1./np.sqrt(2.)))))

        
    return result

def gaussian_norm(x0,s):
    gaussian_norm=0.5
    #gaussian_norm = 0.5 * (1.-np.math.erf((-x0/s)*(1./np.sqrt(2.)))) AIMS
    #gaussian_norm = 0.5 * np.math.erf(x0/(np.sqrt(2)*s)) ME
    return gaussian_norm


# # Calculation parameters

# In[4]:

#Choose colours for each k grid
colours = ("#009E73","#0072B2", "#D55E00", "#CC79A7","#E69F00","#F0E442","#000000","#56B4E9")

#Choose the broadening or range of broadenings to investigate
broads = (np.linspace(0.01,1,int(sys.argv[1])))/hartree

#bin_size of data in nacs_spectrum_kxxx.out
bin_size = 0.01
#maximum energy of data in nacs_spectrum_kxxx.out
#Grid density to discretize the fit on
npoints = 100000
#Number of knots of the piecewise cubic spline fit
nknots = int(sys.argv[2])


# # Parser and calculation

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
print(max_e)
grid_discretization = max_e/npoints
print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
elements = (((dimension*dimension)-dimension)/2)+dimension
print("There are " + str(elements) + " coupling components")
if elements < head_count:
    n_spin = 2
    print("This system is spin unrestricted")
tensor_fit = np.zeros((len(broads),dimension,dimension),dtype=complex)
tensor_norm = np.zeros((len(broads),dimension,dimension),dtype=complex)
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
x =np.linspace(0,max_e,npoints)/hartree
    
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
            
            model_6 = get_natural_cubic_spline_model(e_axis, real, minval=min(e_axis), maxval=max(e_axis), n_knots=nknots)
                
            #NORMAL       
            s_counter = -1
            for s in broads:
                s_counter+=1
                #print('Broadening =',s)
                friction=evaluate_delta_function_2(e_axis, real+1j*imag, 0, s, 'gaussian')*bin_size/hartree
                #friction=evaluate_delta_function(e_axis, real, 0, s, 'gaussian')
                tensor_norm[s_counter,i,j]=friction
                tensor_norm[s_counter,j,i]=np.conj(friction)   
                
            zpf[i,j] = model_6.predict(x)[0]
            zpf[j,i] = model_6.predict(x)[0]
                        
            #FIT
            s_counter = -1
            for s in broads:
                s_counter+=1
                friction=evaluate_delta_function_2(x, model_6.predict(x), 0, s, 'gaussian')*grid_discretization/hartree
                tensor_fit[s_counter,i,j]=friction
                tensor_fit[s_counter,j,i]=friction


# # Plot

# In[10]:


lifetimes_norm = np.zeros((len(broads),dimension))
relaxs_norm = np.zeros((len(broads),dimension))

lifetimes_fit = np.zeros((len(broads),dimension))
relaxs_fit = np.zeros((len(broads),dimension))

zpf_relaxs = np.zeros((dimension))
zpf_lifetimes = np.zeros((dimension))

fig4, ax = plt.subplots(dimension, 1, sharex='all')#, sharey='all')
    
s_counter = -1
for s in broads:
    s_counter +=1
    relaxs_fit[s_counter,:],eigvecs=np.linalg.eigh(tensor_fit[s_counter,:,:])
    #lifetimes_fit[grid_counter,s_counter,:] = np.sort(ps/relaxs_fit[grid_counter,s_counter,:])
    lifetimes_fit[s_counter,:] = ps/relaxs_fit[s_counter,:]
        
    relaxs_norm[s_counter,:],eigvecs=np.linalg.eigh(tensor_norm[s_counter,:,:])
    #lifetimes_norm[grid_counter,s_counter,:] = np.sort(ps/relaxs_norm[grid_counter,s_counter,:])
    lifetimes_norm[s_counter,:] = ps/relaxs_norm[s_counter,:]
    #for i in range(dimension):
    
zpf_relaxs[:],eigvecs=np.linalg.eigh(zpf[:,:])
zpf_lifetimes[:] = ps/zpf_relaxs[:]

for i in range(dimension):
    ax[i].plot(broads*hartree,lifetimes_norm[:,i],'-',linewidth = 1.5,color = colours[0],
            label='GB scheme')
        
    ax[i].plot(broads*hartree,lifetimes_fit[:,i],'--',linewidth = 1.5,color = colours[0],
            label='Interpolation scheme')
    
    ax[i].plot(0,zpf_lifetimes[i],'.',markersize=10,color = colours[0],
            label='Interpolation with no GB')
    ax[i].set_xlim(0.0,broads[-1]*hartree)
    
    if i == 0:
        ax[i].legend(fontsize=10,fancybox=True,framealpha=0,loc=0,ncol=2)
    ax[i].tick_params(labelsize=12)

fig4.set_figheight(12)
fig4.set_figwidth(7)
fig4.text(0.5, 0.00, "Broadening width / eV", ha='center',fontsize=15)
fig4.text(0.02, 0.5, r'Friction mode lifetime / $ps$', va='center', rotation='vertical',fontsize=15)
fig4.savefig('sigma_convergence_nknots'+str(nknots)+'.pdf',transparent=True,bbox_inches='tight')
fig4.savefig('sigma_convergence'+str(nknots)+'.svg',transparent=True,bbox_inches='tight')


# # Frequency-dependence
