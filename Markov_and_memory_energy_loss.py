#!/usr/bin/env python
# coding: utf-8

# # Preamble

# In[1]:


from ase.io import read,write,Trajectory
from ase.visualize import view
import numpy as np
from ase import Atoms
import numpy as np
import os
from ase.db import connect
# coding: utf-8
from ase.calculators.aims import Aims
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from scipy import integrate
from ase.io.trajectory import TrajectoryWriter
import sys
from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import math
from scipy.interpolate import interp1d
from scipy.ndimage import convolve1d


# In[2]:


from ase.units import _hbar, J, s, fs
hbar = _hbar *J
print(hbar) #ev S


# In[3]:


def string2array(string):
    """
    
    """
    dimension = len(string.split(']\n'))
    onedarray = np.fromstring((string.replace('[',' ').replace(']\n',' ')),dtype=float,sep=' ')
    return onedarray.reshape(dimension,dimension)

def shortest_displacement(previous_atoms,next_atoms):
    """
    Calculates the shortest displacement vector between two atoms objects. My solution (albeit flawed)
    to the crossing the unit cell boundary problem. Finds the difference between position vector for each atom
    and the difference between the respectivr lattice vector and difference1. Uses whatever is smallest. 
    Returns that displacement vector
    """
    
    assert(len(previous_atoms)==len(next_atoms))
    unit_cell_lengths = atoms.get_cell_lengths_and_angles()[0:3]
    
    displacement_vector = np.zeros((len(previous_atoms),3))
    
    for i in range(len(previous_atoms)):
        a = next_atoms.positions[i]-previous_atoms.positions[i] 
        b = unit_cell_lengths - abs(a)
        
        for ii in range(3):
            if abs(a[ii]) <= abs(b[ii]):
                displacement_vector[i,ii]=a[ii]
            else:
                displacement_vector[i,ii]=b[ii]
    return(displacement_vector)

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
#    print("Friction max energy = "+str(max_e))
#    print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
    elements = (((dimension*dimension)-dimension)/2)+dimension
#    print("There are " + str(elements) + " coupling components")
    if elements < head_count:
        n_spin = 2
#        print("This system is spin unrestricted")

    bins=np.zeros((int(max_e/discretization)+1))
#    print(len(bins))
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


def inter_bin(t,bins,ys):
        
    if t in bins:
        i = np.where(bins == t)
        return(ys[i][0])
    else:    
        for z,b in enumerate(bins):
            if b<(t):
                lower_bin=b                
                y1 = ys[z]
            if b>(t):
                upper_bin=b        
                y2 = ys[z]
                break
        m = (y2-y1)/(upper_bin-lower_bin)
        c = y1-(m*lower_bin)
        y = m*t+c
        return(y)
    
def man_convolution(x,func1,func2):
    
    """
    Carries out manual convoluton of two functions assumes both 
    are on the same x grid and the grid is uniformly spaced. 
    Normalises with the spacing between each point. Allowed to be complex.
    """
    
    if len(func1)!=len(func2):
        raise ValueError('functions must be same length')
        
    if len(func1)!=len(x):
        raise ValueError('x grid must correspond to functions')
        
    y = np.zeros_like(func1)
    
    for i,t in enumerate(x):
        for ii,tau in  enumerate(x):
        #    
            if tau>t:
                continue
                
            y[i] += func1[ii]*func2[i-ii]
            
        #y[i] = np.sum(func1[0:i]*np.flip(func2[0:i]))
                
    return y*(x[1]-x[0])


# # Input

# In[4]:


con = connect('database.db')
friction_indices = np.array([16,17,18,19]) #counting from 0
final = 100
dimension = 12
n_tensors = 1
## Non-Markovian parameters

cutoffs = np.linspace(0.1,2.5,8)

atoms=con.get_atoms(id=1)


# # Markovian (no memory)

# In[5]:


W = np.zeros((n_tensors,final))
markov_forces = np.zeros((final,len(atoms),3))
all_velocities = np.zeros((final,len(friction_indices),3))
all_tensors = np.zeros((n_tensors,final,dimension,dimension))
for ft in range(1):
    for ts in range(1,final+1):
        atoms = con.get_atoms(id=ts)
        #skip systems that for some reason didnt have pbc:
        if atoms.pbc[0] == False:
            continue
            #atoms.set_cell(cell)
            #atoms.set_pbc(pbc)
        
        #if ts in [474,475,476,477,478,479,480,766,767]:
        #    continue
        
        row = con.get(id=ts)
        try:
            tensor = string2array((row.get('ft'+str(ft+1))))
        except: 
             tensor = np.zeros((dimension,dimension))
         
        # skip finite difference bugged systems
        if np.amax(tensor) > 10:
            continue
            

        vel = atoms.get_velocities() #velocities in sqrt(eV / amu)
        masses = atoms.get_masses()
        force_vec = np.zeros_like(vel)
        if ft ==0:
            all_velocities[ts-1,:,:] = vel[friction_indices,:]
        all_tensors[ft,ts-1,:,:]=tensor
        
        tensor = tensor * 1e12 # s-1
        tensor = tensor * (1/s) #inverse(ASE_TIME)
        
        ##################################
        i_cart = -1
        i_atom = 0
        for i in range(dimension):
            i_cart = i_cart + 1
            if (i_cart>2):
                i_cart = 0
                i_atom = i_atom + 1
        #            
            j_atom = 0
            j_cart = -1          
            for j in range(dimension):
                j_cart = j_cart + 1
                if (j_cart>2):
                    j_cart = 0
                    j_atom = j_atom + 1
                mass_factor=np.sqrt(masses[friction_indices[i_atom]])*np.sqrt(masses[friction_indices[i_atom]])
                tensor[i,j]=tensor[i,j]*mass_factor
                #force_vec[friction_indices[i_atom],i_cart]+= tensor[i,j]*vel[friction_indices[j_atom],j_cart]*mass_factor         
    
        tensor_full = np.zeros((len(vel)*3,len(vel)*3))
        x = (3*len(vel))-tensor.shape[0]
        tensor_full[x:x+tensor.shape[0], x:x+tensor.shape[1]] = tensor
        force_vec=np.dot(tensor_full,vel.flatten())


        W[ft,ts-1] += np.dot(vel.flatten(),np.dot(tensor_full,vel.flatten()))*2*fs
        #W[ft,ts-1] += np.dot(displacement_vec.flatten(),force_vec)
        #for ii in range(66):
        markov_forces[ts-1,:,:] = force_vec.reshape(len(atoms),3)
            #W[ft,ts-1] += np.dot(displacement_vec[ii],force_vec[ii])
    


# ## Plots

# In[6]:


plt.figure(figsize=(10,6))
atom_labels=['Cx','Cy','Cz','Ox','Oy','Oz']
lines = ['-','-','-','--','--']
time_axis = np.arange(0,final*2,2)
for ft in range(n_tensors):
    #if ft in [0,1,2,5]:
    #    continue
    for i in range(6):
        for j in range(6):
            if j<i:
                continue
            if np.average(np.abs(all_tensors[0,:,i,j]))>0.1:
                print(str(ft)+atom_labels[i]+atom_labels[j])
                
                plt.plot(time_axis,all_tensors[ft,:,i,j],label=atom_labels[i]+atom_labels[j],
                        linestyle=lines[ft])

plt.legend()
#plt.xlim(0,200)#
plt.xlabel('Time / fs',fontsize=20)
plt.ylabel(r'$\Lambda_{ij}$ / $ps^{-1}$',fontsize=20)


# In[7]:


plt.figure(figsize=(10,6))
atom_labels=['C','O']
cart_labels = ['x','y','z']
lines = ['-','--']
for atom in range(2):
    for cart in range(3):

        plt.plot(time_axis,all_velocities[:,atom,cart],label=str(atom_labels[atom])+str(cart_labels[cart]),
                linestyle=lines[atom])

plt.legend()
#plt.xlim(0,600)
plt.xlabel('Time / fs',fontsize=20)
plt.ylabel(r'Velocity / $\sqrt{\frac{eV}{amu}}$',fontsize=20)


# In[8]:


plt.figure(figsize=(10,10))
colours=['red','blue','green','purple','grey','cyan']
labels = ['\sigma = $0.6 eV','\sigma=$0.1 eV','\sigma=$0.01 eV','\ \Lambda(E_k)$ linear','\ \Lambda(E_k)$ cubic','\ $mem-MDEF']
lines=['-','-','-','--','--','-.']
for ft in range(n_tensors):
    plt.plot(time_axis,W[ft,:],label=r'$'+labels[ft],color=colours[ft],linestyle=lines[ft])
    #plt.plot(time_axis,integrate.cumtrapz(W[ft,:],time_axis,initial=0))
    plt.hlines(np.average(W[ft,:]),0,np.max(time_axis),color=colours[ft],label='average')
    plt.legend()
    #plt.xlim(0,200)
    #plt.ylim(-0.0001,0.0005)
plt.xlabel('Time / fs',fontsize=20)
plt.ylabel('Work / eV',fontsize=20)


# # Non-Markovian
# 

# For each cutoff:
#     1. Need to read in friction_memory_kernel.out per time
#     2. Fourier transform
#     3. Interpolate
#     4. Convolute to get eta_t per ts
#     

# In[9]:


## Calculate derived parameters that are same for every step (e.g frequencies & times & mass factors)
frequency_list = []
times_list = []
eta_bar_t_list = []

for co,cutoff in enumerate(cutoffs):
    ts=1
    row = con.get(id=ts)
    path = str(sys.argv[1])+str(ts)+'/friction_memory_kernel.out'
    bins,re_memory_kernel,im_memory_kernel,dimension,max_e = read_memory_kernel(path)
    masses = atoms.get_masses() # amu
    
    frequencies=((bins[bins<cutoff])/hbar)/s #ase units inverse time
    times = 1/frequencies #ase units time
  
    times = np.append(times,0)

            
    frequency_list.append(frequencies)
    times_list.append(times)
    eta_bar_t = np.zeros((final,dimension,dimension,len(times)))
    eta_bar_t_list.append(eta_bar_t)
    
    


# ## PARSE ALL DATA

# In[10]:


if 'read' in sys.argv[1:]:
    print('reading')
    raw_data = np.load('raw_memory.npy')

else:
    raw_data = np.zeros((final,dimension,dimension,len(bins)))
    for ts in range(1,final+1):
        row = con.get(id=ts)
        try:
            atoms = con.get_atoms(id=ts)
        except:
            print('cannot get atoms for '+str(ts)+' - continuing')
            continue
    
        path = str(sys.argv[1])+str(ts)+'/friction_memory_kernel.out'
        try:
            bins,re_memory_kernel,im_memory_kernel,dimension,max_e = read_memory_kernel(path)
        except:
            print('cannot get mem_kernel for '+str(ts)+' - continuing')
            continue
                  
        raw_data[ts-1,:,:,:] = re_memory_kernel
    np.save('raw_memory.npy',raw_data)
# ## Fourier transform

# In[11]:


for ts in range(1,final+1):
#for ts in range(1,101):
    print(ts)
    row = con.get(id=ts)
    try:
        atoms = con.get_atoms(id=ts)
    except:
        print('cannot get atoms for '+str(ts)+' - continuing')
        continue
    
    
    re_memory_kernel=raw_data[ts-1,:,:,:]
    
    for co,cutoff in enumerate(cutoffs):
        times = times_list[co]
        frequencies = frequency_list[co]
        eta_bar_t = eta_bar_t_list[co]

            
        for i in range(dimension):
            i_atom = int(i/3)
            for j in range(dimension):
                if j<i:
                    continue
                
                
                j_atom = int(j/3)
                mass_factor=np.sqrt(masses[friction_indices[i_atom]])*np.sqrt(masses[friction_indices[i_atom]]) #amu
                lambda_omega = re_memory_kernel[i,j,:]*mass_factor/(fs*1000)
        
        
                func = np.zeros(len(frequencies))

                for t in range(len(times)):

                    #if times[t]==0:
                    #    eta_bar_t[ts-1,i,j,t]=np.average(lambda_omega)
                    #    continue

                    #for w in range(len(frequencies)):

                        #func[w]=(lambda_omega[w])*np.cos(frequencies[w]*times[t])
                        

                    func=(lambda_omega[0:len(frequencies)])*np.cos(frequencies*times[t])
                    func[0]=0
                # amu ase(time)^-2
                    eta_bar_t[ts-1,i,j,t]=np.trapz(func[:],frequencies[:])
        eta_bar_t_list[co]=eta_bar_t
    #eta_bar_t_list.append(eta_bar_t)
        


# ## Interpolation

# In[12]:


times_up_list = []
times_sort_list = []
eta_bar_inter_list = []
eta_bar_sort_list = []

#times_up = np.linspace(-10*fs,10*fs,500)
times_up = np.linspace(-9*fs,9*fs,500)
below_zero = np.argwhere(times_up < 0.0)
for co,cutoff in enumerate(cutoffs):
    
    #get info
    times = times_list[co]
    frequencies = frequency_list[co]    
    eta_bar_t = eta_bar_t_list[co]
    
    #sort times
    time_sidx = np.argsort(times)

    times_sort = (times)[time_sidx]
    times_sort_list.append(times_sort)
    
    
    times_up_list.append(times_up)
    
    eta_bar_inter = np.zeros((final,dimension,dimension,len(times_up)),dtype=complex)
    eta_bar_sort = np.zeros((final,dimension,dimension,len(times_sort)))
    
    
    
    for ts in range(1,final+1):
    #for ts in range(1,101):
        for i in range(dimension):       
            for j in range(dimension):
                if j<i:
                    continue  
            
                unsort = eta_bar_t[ts-1,i,j,:]
            
                eta_bar_sort[ts-1,i,j,:]=unsort[time_sidx]
            
                f = interpolate.interp1d(times_sort,unsort[time_sidx],fill_value="extrapolate")
                
                eta_bar_inter[ts-1,i,j,:] = f(times_up)
                
                eta_bar_inter[ts-1,i,j,below_zero] = 0
    
    eta_bar_sort_list.append(eta_bar_sort)    
    eta_bar_inter_list.append(eta_bar_inter)           


# ## Convolution
# 

# In[13]:


eta_t_list = []
bugged_ones = [275,397,649,766,767,786]
for co,cutoff in enumerate(cutoffs):
    print(co)
    #get info
    times_up=times_up_list[co]    
    eta_bar_inter = eta_bar_inter_list[co]
    dx = times_up[1]-times_up[0]
    #find cutoff

    cutoff_freq = (cutoff/hbar)/s
    #factors
    sinc = cutoff_freq*np.sinc((times_up*cutoff_freq)/2)    
    
    #exp_factor = np.cos(0.5*times_up*cutoff_freq) + (1j * np.sin(0.5*times_up*cutoff_freq))
    exp_factor = np.cos(0.5*times_up*cutoff_freq)

    #eta_t = np.zeros((final,dimension,dimension,len(times_up)),dtype=complex)
    eta_t = np.zeros((final,dimension,dimension,len(times_up)),dtype=float)
    
    convolute_factor = sinc*exp_factor
    
    for ts in range(1,final+1):
    #for ts in range(1,101):
        if ts in bugged_ones:
            continue
        for i in range(dimension):    
            for j in range(dimension):
                if j<i:
                    continue  
            

                #eta_t[ts-1,i,j,:] = man_convolution(times_up,eta_bar_inter[ts-1,i,j,:],convolute_factor)
                
                
                eta_t[ts-1,i,j,:] = np.convolve(eta_bar_inter[ts-1,i,j,:],convolute_factor,'same')*dx
                
                #eta_t[ts-1,i,j,:] = convolve1d(np.real(eta_bar_inter[ts-1,i,j,:]),convolute_factor)*dx
                
                #eta_t[ts-1,i,j,:] = signal.convolve(eta_bar_inter[ts-1,i,j,:],convolute_factor,'same')*dx

            
            
    
    eta_t_list.append(eta_t)



# In[14]:


plt.plot(times_up,convolute_factor)
plt.plot(times_up,eta_bar_inter[0,0,0,:])
plt.plot(times_up,eta_t[0,0,0,:])


# ## x Velocity to get force and x displacement vector to get work then integrate

# In[15]:


## Interpolate velocities for all possible t_primes

## Fairly bad at the moment - Can make this a lot more efficient by actually interpolating the entire data range
## i.e to some spline fit in  one go
#denser_ts = np.linspace(1,final,(final*10)+1,endpoint=True)
denser_ts = np.arange(1,final,0.05)
denser_ts = np.append(denser_ts,final)
denser_ts = np.round(denser_ts,2)
print(denser_ts)
velocities_inter = np.zeros((len(denser_ts),int(dimension/3),3))


for idx,dts in enumerate(denser_ts):
    
    if dts in range(1,final+1):
         velocities_inter[idx,:,:] = all_velocities[int(dts-1),:,:]
    
    else:
        
        low_ts = int(math.floor(dts))
        high_ts = int(math.ceil(dts))
        
#         print(low_ts)
#         print(dts)
#         print(high_ts)
#         print('-----')
        
        if high_ts ==final+1:
            high_ts=final
        
        fst = (all_velocities[low_ts-1,:,:]).flatten()
        snd = (all_velocities[high_ts-1,:,:]).flatten()
        
        linfit = interp1d([0,1], np.vstack([fst, snd]), axis=0) # unnecsairly repeated for all points between
        # two bins
        
        velocities_inter[idx,:,:]=(linfit(dts-low_ts)).reshape(int(dimension/3),3)        


# ### Plot check interpolated velocities

# In[16]:


plt.figure(figsize=(10,6))
atom_labels=['C','O']
cart_labels = ['x','y','z']
lines = ['-','--']
time_axis = np.arange(0,final*2,2)
for atom in range(2):
    for cart in range(3):

        plt.plot((denser_ts-1)*2,velocities_inter[:,atom,cart],label=str(atom_labels[atom])+str(cart_labels[cart]),
                linestyle=lines[atom])
        
        plt.plot(time_axis,all_velocities[:,atom,cart],'.',label=str(atom_labels[atom])+str(cart_labels[cart]))

plt.legend()

plt.xlim(0,100)
plt.xlabel('Time / fs',fontsize=20)
plt.ylabel(r'Velocity / $\sqrt{\frac{eV}{amu}}$',fontsize=20)


# ### Mem integral

# In[68]:


mem_cutoff = 7.5 #fs
max_t = mem_cutoff*fs #upper limit for integration
run = True
nm_work=np.zeros((len(cutoffs),final))
force_vec_all = np.zeros((len(cutoffs),final,len(atoms),3))
for co,cutoff in enumerate(cutoffs):
    
    if cutoff<0.5:
        continue
    print(co)
    eta_t = eta_t_list[co]
    times=times_up_list[co]
    
    ##################################t s loop #####################################################
    for ts in range(1,final+1):
    #for ts in range(1,500+1):
        atoms = con.get_atoms(id=ts)
        if atoms.pbc[0] == False:
            continue
        if ts in bugged_ones:
            continue

        
        # WANT t primes as denser ts up to and including ts
  
        t_primes = ((denser_ts[denser_ts<=ts])-1)*2

        #Integrate dt' eta(t-t')vel(t') between t0 and t
        time_step = (ts-1)*2*fs #ts in ase time
        time_axis = time_step - (t_primes*fs)
        integrand = np.zeros([len(t_primes),dimension,dimension])
        #############################################t prime loop########################################  
        for tp,t_prime in enumerate(t_primes):  #tp is the index of denser time, t' where evaluating integral
            # t_primes do not neccesarily correspond to actual time steps calculated. 
            if ts ==1:
                integrand[0,:,:]=0
                continue

            real_below_ts = int((math.floor(t_prime)/2)+1)            
            if real_below_ts in bugged_ones:
                continue
            
            
            #Need a way to know which t_primes correspond to actual calcs
            #And the others are just interpolated between the nearest two  real calcs        
            #print('val1 ' + str(((ts-1)*2)))

            if time_step-(t_prime*fs) > max_t:
                integrand[tp,:,:]=0
                continue
            if t_prime == ((ts-1)*2):
                integrand[tp,:,:]=eta_t[real_below_ts-1,:,:,0]
                continue  
            
            

            i_cart = -1
            i_atom = 0
            for i in range(dimension):
                i_cart = i_cart + 1
                if (i_cart>2):
                    i_cart = 0
                    i_atom = i_atom + 1
            #            
                j_atom = 0
                j_cart = -1          
                for j in range(dimension):
                    j_cart = j_cart + 1
                    if (j_cart>2):
                        j_cart = 0
                        j_atom = j_atom + 1
                    
                    if j<i:
                        continue


                    
                    eta_t_sort = eta_t[real_below_ts-1,i,j,:]*-1 #ASSUMPTION: eta_t function CHANGES SLOWLY
                    
                    #Should really interpolate eta_t[real_below] and eta_[real_above] to time in between
                    
                    
                    
                    #dont really need tis interpolation if ets t is on a uniform grid? make
                    #t_primes have the same spcing as times_up
                    integrand[tp,i,j] = inter_bin(time_step-(t_prime*fs),times,eta_t_sort)
                    #integrand[tp,i,j]*=vel[friction_indices[j_atom],j_cart]
                    integrand[tp,i,j]*=velocities_inter[tp,j_atom,j_cart]
                    integrand[tp,j,i]= integrand[tp,i,j]

        ########### END t prime loop  
        vel = atoms.get_velocities()
        
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

                force_vec_all[co,ts-1,friction_indices[i_atom],i_cart]+=np.trapz(integrand[:,i,j],time_axis)

            nm_work[co,ts-1]+=np.dot(vel[friction_indices[i_atom],i_cart],force_vec_all[co,ts-1,friction_indices[i_atom],i_cart])

#         if ts == 1:
#             ts_previous = 1
#             ts_next = ts+1
#         elif ts == final:
#             ts_next = final
#             ts_previous = ts-1
#         else:
#             ts_previous = ts-1
#             ts_next = ts+1
        
        
        
        #displacement_vec = shortest_displacement(con.get_atoms(id=ts_previous),con.get_atoms(id=ts_next))/2
        #dot product force and displacement now units is eV = energy
        #for ii in range(len(atoms)):
        #    nm_work[co,ts-1] += np.dot(displacement_vec[ii,:],force_vec_all[co,ts-1,ii,:])
        #for ii in range(len(atoms)):
         #   nm_work[co,ts-1] += np.dot(force_vec_all[co,ts-1,ii,:],vel[ii,:])*2*fs
nm_work = nm_work*2*fs


# # Plot

# In[89]:


colours = ['firebrick','red','orange','skyblue','blue','purple','grey','green','yellow','pink','black','fuchsia','gold']

labels = ['\sigma = $0.6 eV','\sigma=$0.1 eV','\sigma=$0.01 eV','\ \Lambda(\mathrm{E_k})$ linear','\ \Lambda(\mathrm{E_k})$ cubic','\int\eta({\\tau})v(t\')\ dt\'\ $linear']
fig3, ax = plt.subplots(1, 1, sharex='all', sharey='all'
                      )
time_axis = np.arange(0,final*2,2)
for ft in range(n_tensors):
    #ax.plot(time_axis,np.cumsum(np.sum(B[ft,:,:],axis=1)),'o',markersize = 2,mfc='none',label=labels[ft])
    ax.plot(time_axis,integrate.cumtrapz(W[ft,:],time_axis,initial=0),'--',
            color=colours[ft],
            markersize = 2,mfc='none',label=r'$'+labels[ft])
    #y_position = integrate.cumtrapz(np.sum(B[ft,:,:],axis=1),time_axis,initial=0)[-1]
    #ax.text(time_axis[-125],y_position, r'$'+labels[ft],
    #        rotation=5,
    #       fontsize=12)

for co,cutoff in enumerate(cutoffs):
    #if cutoff<1:
    #    continue
    ax.plot(time_axis,integrate.cumtrapz(nm_work[co,:],time_axis,initial=0),'-',
            color=colours[co+1],markersize = 2,mfc='none',label=r'$\omega_c = {:.2f}$ eV'.format(cutoff))
    
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlim(0,final*2)
#ax.set_xlim(0,500*2)
#ax.set_ylim(0,0.02)
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.legend(fontsize=15,fancybox=True,framealpha=0,loc=0,ncol=1)
fig3.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig3.text(-0.02, 0.5, r'Cumulative energy loss / eV', va='center', rotation='vertical',fontsize=15)
fig3.set_figheight(10)
fig3.set_figwidth(7)
fig3.savefig('e_loss_mem_cutoff.pdf',transparent=True,bbox_inches='tight')


# In[71]:


fig5, ax = plt.subplots(1, 1, sharex='all', sharey='all')

for ft in range(n_tensors):
    #ax.plot(time_axis,np.cumsum(np.sum(B[ft,:,:],axis=1)),'o',markersize = 2,mfc='none',label=labels[ft])
    ax.plot(time_axis,W[ft,:],'--',
            color=colours[ft],
            markersize = 2,mfc='none',label=r'$'+labels[ft])


for co,cutoff in enumerate(cutoffs):
    #if cutoff<1:
    #    continue
    ax.plot(time_axis,nm_work[co,:],'-',
            color=colours[co+1],markersize = 2,mfc='none',label=r'$\omega_c = {:.2f}$ eV'.format(cutoff))
    
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlim(0,200)
#ax.set_ylim(0,0.02)
#ax.xaxis.set_minor_locator(MultipleLocator(50))
#ax.yaxis.set_minor_locator(MultipleLocator(0.01))
ax.legend(fontsize=15,fancybox=True,framealpha=0,loc=0,ncol=1)
fig5.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig5.text(-0.02, 0.5, r'Energy loss / eV', va='center', rotation='vertical',fontsize=15)
fig5.set_figheight(10)
fig5.set_figwidth(7)
fig5.savefig('work_mem_cutoff.pdf',transparent=True,bbox_inches='tight')


# In[72]:


fig4, ax = plt.subplots(1, 1, sharex='all', sharey='all')

atom = 2
for co,cutoff in enumerate(cutoffs):
    ax.plot(time_axis,force_vec_all[co,:,-atom,2],linewidth=1,color=colours[co+1])

ax.plot(time_axis,markov_forces[:,-atom,2], '--',linewidth=3,color=colours[0])
ax.set_xlim(0,200)
fig4.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig4.text(-0.02, 0.5, r'Force/ ', va='center', rotation='vertical',fontsize=15)
fig4.set_figheight(7)
fig4.set_figwidth(7)                      
fig4.savefig('forces.pdf',transparent=True,bbox_inches='tight')


# In[73]:


print(atoms[-2])


# In[80]:


for co,cutoff in enumerate(cutoffs):
    y_data = eta_bar_t_list[co]
    x_data = times_list[co]
    for ts in range(1,2):
        plt.plot(x_data/fs,y_data[ts-1,2,2,:],color=colours[co+1],marker='.')
        plt.xlim(0,mem_cutoff+(2))
plt.vlines(mem_cutoff,-1,2)


# In[81]:


for co,cutoff in enumerate(cutoffs):
    y_data = eta_bar_sort_list[co]
    x_data = times_sort_list[co]
    for ts in range(1,2):
        plt.plot(x_data/fs,y_data[ts-1,2,2,:],color=colours[co+1],marker='.')
        plt.xlim(0,mem_cutoff+(2))
plt.vlines(mem_cutoff,-1,2)
plt.savefig('eta_bar.pdf')


# In[82]:


for co,cutoff in enumerate(cutoffs):
    y_data = eta_bar_inter_list[co]
    x_data = times_up_list[co]
    for ts in range(1,2):
        plt.plot(x_data/fs,y_data[ts-1,2,2,:],color=colours[co+1],marker='.')
        plt.xlim(0,mem_cutoff+(2))
plt.vlines(mem_cutoff,-1,2)
plt.savefig('eta_bar_inter.pdf')


# In[83]:


for co,cutoff in enumerate(cutoffs):
    y_data = eta_t_list[co]
    x_data = times_up_list[co]
    for ts in range(1,final+1):
        plt.plot(x_data/fs,y_data[ts-1,2,2,:],color=colours[co+1])
        
plt.vlines(0,-5,5)
plt.vlines(mem_cutoff,-5,5)
plt.xlim(0,mem_cutoff+(2))
plt.savefig('eta_t.png',dpi=300)


# In[ ]:




