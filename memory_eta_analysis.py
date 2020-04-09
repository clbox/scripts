import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
from ase.io import read
import math
hbar = _hbar * J * s 
ps = fs*1000

def read_memory_kernel(path,treat_complex=True):
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

    elements = int((((dimension*dimension)-dimension)/2)+dimension)
    if elements < head_count:
        n_spin = 2 
        #print("This system is spin unrestricted")

    bins=np.zeros((int(max_e/discretization)+1))
    re_memory_kernel=np.zeros((elements,len(bins)))
    im_memory_kernel=np.zeros_like(re_memory_kernel)
    e=0
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                c=0
                e +=1
            if any(x in line for x in header):
                continue
            else:
                re_memory_kernel[e-1,c]=float(line.split()[1])
                if treat_complex:
                    im_memory_kernel[e-1,c]=float(line.split()[2])
                bins[c]=float(line.split()[0])
                c +=1

    return(bins,re_memory_kernel,im_memory_kernel,dimension,max_e)

def frequency_interpolate(bins,re_memory_kernel,dimension,mem_cutoff):
    """Add extra bins in frequency domain"""
    
    time_bins = np.linspace(4*fs,mem_cutoff,500)
    add_bins = (1/(time_bins))*hbar #eV
    new_bins = np.append(bins,add_bins)
    new_bins = np.sort(new_bins)

    elements = np.shape(re_memory_kernel)[0]
    new_data = np.zeros((elements,len(new_bins)))


    for e in range(elements):
        f = interp1d(bins,re_memory_kernel[e,:])
        new_data[e,:] = f(new_bins)

    return new_bins,new_data

def generate_domains(new_bins):
    frequencies=new_bins/hbar #ase units inverse time
    times = 1/frequencies #ase units time

    times = np.append(times,0)
    times = np.sort(times) #check consequencies
    times[-1] = 2**32 #check consequencies assumes 0 freq

    return times,frequencies

def get_friction_masses(path,friction_indices):
    """ get friction masses from geometry.in"""


    atoms = read(path)

    masses = atoms.get_masses()

    friction_masses = masses[friction_indices]

    return friction_masses

def element_index(dimension):
    "Returns list of indices for flattened triangular array"
    elements = int((((dimension*dimension)-dimension)/2)+dimension)
    a = np.arange(0,elements)
    i=0 
    j=0
    indices = []
    for e in a:
        indices.append((i,j))
        if j == dimension-1: 
            i += 1 
            j = i 
            continue
        j+=1
    return(indices)

def fourier_transform(times,masses,frequencies,new_data,dimension):
    """Fourier transform to time domain"""
    elements = np.shape(new_data)[0]
    indices = element_index(dimension)

    eta_bar_t = np.zeros((elements,len(times)))
    func = np.zeros(len(frequencies))
    cos_factor = np.cos(frequencies*times[:,None])


    lambda_omega = new_data[:,0:len(frequencies)] #convert from ps-1
        
    for e in range(elements):
        i_atom,j_atom = (indices[e])[0] // 3, (indices[e])[1] // 3
        func = lambda_omega[e,None,:] * cos_factor * np.sqrt(masses[i_atom]*masses[j_atom])
        func[:,0]=0
        eta_bar_t[e,:]=np.trapz(func,frequencies,1)


    return eta_bar_t

def time_interpolate(times,eta_bar_t,time_step,mem_cutoff):
    """At the moment it is neccesary to interpolate in the  time domain
    as well to get a uniform distribution of bins for the convolution 
    
    Interpolated time domain has limits derived from both the memory cutoff
    and the time step. the spacing of the interpolated domain is chosen to ensure
    it can match the interpolated 'nuclear time domain'. i.e the spacing is a
    factor of the max nuclear time.
    
    """

    elements = np.shape(eta_bar_t)[0]

    dt = (time_step/40)
    n_points = math.ceil(mem_cutoff/dt)
    time_e_max = dt * (n_points-1)
    n_points = n_points + n_points-1
    times_forward = np.linspace(0,time_e_max,n_points)
    times_up = np.append(times_forward,-times_forward[1:])
    times_up = np.sort(times_up)
    below_zero = np.argwhere(times_up < 0.0)

    eta_bar_inter = np.zeros((elements,len(times_up)))         

    for e in range(elements):
        f = interp1d(times,eta_bar_t[e,:],fill_value="extrapolate")
        
        eta_bar_inter[e,:] = f(times_up)
        
        eta_bar_inter[e,below_zero] = 0
    

    return times_up,eta_bar_inter    

def convolute(times_up,eta_bar_inter,cutoff):
    """Convolute with sinc factor in time domain"""

    elements = np.shape(eta_bar_inter)[0]
    dx = times_up[1]-times_up[0]
    #find cutoff

    cutoff_freq = (cutoff/hbar)
    #factors
    sinc = cutoff_freq*np.sinc((times_up*cutoff_freq)/2)    
    
    #exp_factor = np.cos(0.5*times_up*cutoff_freq) + (1j * np.sin(0.5*times_up*cutoff_freq))
    exp_factor = np.cos(0.5*times_up*cutoff_freq)

    #eta_t = np.zeros((final,dimension,dimension,len(times_up)),dtype=complex)
    eta_t = np.zeros((elements,len(times_up[times_up >= 0.0])))
    
    convolute_factor = sinc*exp_factor
    
    for e in range(elements):
        eta_t[e,:] = (np.convolve(eta_bar_inter[e,:],convolute_factor,'same')*dx)[times_up >= 0.0]
    return eta_t

def get_eta_t(path,mem_cutoff,time_step,masses):

    mem_cutoff = mem_cutoff * fs

    time_step = time_step * fs

    bins,re_memory_kernel,im_memory_kernel,dimension,max_e = read_memory_kernel(path)

    re_memory_kernel /= ps

    new_bins,new_data = frequency_interpolate(bins,re_memory_kernel,dimension,mem_cutoff)

    times,frequencies = generate_domains(new_bins)

    eta_bar_t = fourier_transform(times,masses,frequencies,new_data,dimension)

    times_up,eta_bar_inter = time_interpolate(times,eta_bar_t,time_step,mem_cutoff)

    cutoff = np.max(bins)

    print('Cutoff: '+ str(cutoff) + ' eV')

    eta_t = convolute(times_up,eta_bar_inter,cutoff)

    return times_up[times_up >= 0.0], eta_t, dimension
