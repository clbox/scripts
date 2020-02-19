import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
from ase.db import connect
import math
hbar = _hbar * J


class Postprocessed_memory:

    """ Takes Raw data, array of n_structures,dimension,dimension,energy_bins 
        and calculates the spectral density, Fourier transform, with 
        convolution and calculates energy loss due to friction including
        the effects of memory 
        
        Can do this at multiple cutoff energies if supplied, or one"""
    
    def __init__(self,bins,raw_data,cutoffs,mem_cutoff,friction_indices,time_step,con):      
        self.raw_data = raw_data
        self.dimension = np.shape(raw_data)[1]
        self.cutoffs = cutoffs
        self.mem_cutoff = mem_cutoff * fs
        self.steps = np.shape(raw_data)[0]
        self.friction_indices = friction_indices
        self.time_step = time_step * fs
        self.bins = bins
        self.con = con

    def frequency_interpolate(self):
        """Add extra bins in frequency domain"""
        
        time_bins = np.linspace(4e-15,self.mem_cutoff/fs,500) #fs
        add_bins = (1/(time_bins))*hbar #eV
        new_bins = np.append(self.bins,add_bins)
        new_bins = np.sort(new_bins)

        new_data = np.zeros((self.steps,self.dimension,self.dimension,len(new_bins)))

        for ts in range(self.steps):
                for i in range(self.dimension):
                    for j in range(i,self.dimension):
                        f = interp1d(self.bins,self.raw_data[ts,i,j,:])
                        new_data[ts,i,j,:] = f(new_bins)
        self.new_data = new_data
        self.new_bins = new_bins

    def generate_domains(self):
        frequency_list = []
        times_list = []
        self.eta_bar_t_list = []
        for cutoff in self.cutoffs:
            frequencies=((self.new_bins[self.new_bins<cutoff])/hbar)/s #ase units inverse time
            times = 1/frequencies #ase units time
        
            times = np.append(times,0)
            times = np.sort(times) #check consequencies
            times[-1] = 2**32 #check consequencies assumes 0 freq

            frequency_list.append(frequencies)
            times_list.append(times)
            eta_bar_t = np.zeros((self.steps,self.dimension,self.dimension,len(times)))
            self.eta_bar_t_list.append(eta_bar_t)

        self.frequency_list = frequency_list
        self.times_list = times_list
    
    def fourier_transform(self):
        """Fourier transform to time domain"""
        friction_indices = self.friction_indices
        dimension =  self.dimension

        try:
            atoms = self.con.get_atoms(id=1)
        except:
            print('cannot get atoms for id = '+str(1))
        masses = atoms.get_masses()[friction_indices]

        for co in range(len(self.cutoffs)):
            times = self.times_list[co]
            frequencies = self.frequency_list[co]
            func = np.zeros(len(frequencies))
            eta_bar_t = self.eta_bar_t_list[co]
            cos_factor = np.cos(frequencies*times[:,None])
            print('max cos: ' + str(np.max(cos_factor)))
            for ts in range(self.steps):
                lambda_omega = self.new_data[ts,:,:,0:len(frequencies)]/(fs*1000) #convert from ps-1
                for i in range(dimension):
                    i_atom = i // 3       
                    for j in range(i,dimension):
                        j_atom = j // 3
                        lambda_omega*=np.sqrt(masses[i_atom]*masses[j_atom])
                        func = lambda_omega[i,j,None,:] * cos_factor
                        func[:,0]=0
                        print('max func:' + str(np,max(func)))
                        eta_bar_t[ts,i,j,:]=np.trapz(func,frequencies,1)
            print('max etabart:' + str(np,max(eta_bar_t)))
            self.eta_bar_t_list[co]=eta_bar_t

    def time_interpolate(self):
        """At the moment it is neccesary to interpolate in the  time domain
        as well to get a uniform distribution of bins for the convolution 
        
        Interpolated time domain has limits derived from both the memory cutoff
        and the time step. the spacing of the interpolated domain is chosen to ensure
        it can match the interpolated 'nuclear time domain'. i.e the spacing is a
        factor of the max nuclear time.
        
        """
        self.times_up_list = []
        self.eta_bar_inter_list = []


        dt = (self.time_step/10)
        n_points = math.ceil(self.mem_cutoff/dt)
        time_e_max = dt * (n_points-1)
        times_forward = np.linspace(0,time_e_max,n_points)
        times_up = np.append(times_forward,-times_forward[1:])
        times_up = np.sort(times_up)
        below_zero = np.argwhere(times_up < 0.0)

        for co in range(len(self.cutoffs)):
            
            times = self.times_list[co]   
            eta_bar_t = self.eta_bar_t_list[co]
            
            self.times_up_list.append(times_up)
            
            eta_bar_inter = np.zeros((self.steps,self.dimension,self.dimension,len(times_up)))         
            
            for ts in range(self.steps):
                for i in range(self.dimension):
                    for j in range(i,self.dimension):
                    
                        f = interp1d(times,eta_bar_t[ts,i,j,:],fill_value="extrapolate")
                        
                        eta_bar_inter[ts,i,j,:] = f(times_up)
                        
                        eta_bar_inter[ts,i,j,below_zero] = 0
            

            self.eta_bar_inter_list.append(eta_bar_inter)       

    def convolute(self):
        """Convolute with sinc factor in time domain"""
        self.eta_t_list = []

        for co,cutoff in enumerate(self.cutoffs):
            times_up=self.times_up_list[co]    
            eta_bar_inter = self.eta_bar_inter_list[co]
            dx = times_up[1]-times_up[0]
            #find cutoff

            cutoff_freq = (cutoff/hbar)/s
            #factors
            sinc = cutoff_freq*np.sinc((times_up*cutoff_freq)/2)    
            
            #exp_factor = np.cos(0.5*times_up*cutoff_freq) + (1j * np.sin(0.5*times_up*cutoff_freq))
            exp_factor = np.cos(0.5*times_up*cutoff_freq)

            #eta_t = np.zeros((final,dimension,dimension,len(times_up)),dtype=complex)
            eta_t = np.zeros((self.steps,self.dimension,self.dimension,len(times_up)))
            
            convolute_factor = sinc*exp_factor
            
            for ts in range(self.steps):
                for i in range(self.dimension):
                    for j in range(i,self.dimension):
                        eta_t[ts,i,j,:] = np.convolve(eta_bar_inter[ts,i,j,:],convolute_factor,'same')*dx 
            self.eta_t_list.append(eta_t)

    def get_velocities(self):
        
        self.all_velocities = np.zeros((self.steps,len(self.friction_indices),3))
        for i in range(self.steps):
            atoms = self.con.get_atoms(id=i+1)
            self.all_velocities[i,:,:] = atoms.get_velocities()[self.friction_indices,:]

    def velocitiy_interpolation(self):
        """To properly evaluate the memory integral, need to carry out
        integration on a scale denser than the nuclear time step scale"""

        steps = self.steps
        time_step = self.time_step
        dt = (self.times_up_list[0])[1]-(self.times_up_list[0])[0]
        max_t = (steps-1)*time_step
        n_old = steps
        n_new = int(round(max_t/dt))+1
        inter_time_scale = np.linspace(0,max_t,n_new)
        old_time_scale = np.linspace(0,max_t,n_old)
        self.velocities_inter = np.zeros((len(inter_time_scale),len(self.friction_indices),3))

        for atom in range(len(self.friction_indices)):
            for cart in range(3):
                linfit = interp1d(old_time_scale,self.all_velocities[:,atom,cart]) 

                self.velocities_inter[:,atom,cart] = linfit(inter_time_scale)

        self.inter_time_scale = inter_time_scale

    def mem_integral(self):   
        old_time_scale = np.arange(0,self.steps*self.time_step,self.time_step)
        dimension = self.dimension
        mem_cutoff = self.mem_cutoff
        inter_time_scale = self.inter_time_scale
        velocities_inter = self.velocities_inter
        force_vec = np.zeros((len(self.cutoffs),len(inter_time_scale),len(self.friction_indices),3))
        nm_work=np.zeros((len(self.cutoffs),len(inter_time_scale)))

        for co in range(len(self.cutoffs)):

            eta_t = self.eta_t_list[co]*-1
            fit = interp1d(old_time_scale,eta_t[:,:,:,:],kind='linear',axis=0)
            eta_t_fit = fit(inter_time_scale)

            for i in range(dimension):
                i_cart = i % 3
                i_atom = i // 3       
                for j in range(i,dimension):
                    j_cart = j % 3
                    j_atom = j // 3
                    for ts, time_step in enumerate(inter_time_scale):
                        t_primes = inter_time_scale[inter_time_scale<=time_step]
                        time_axis = time_step - t_primes
                        integrand = np.zeros([len(t_primes),dimension,dimension])

                        if ts == 0:
                            integrand[0,i,j]=0
                            continue
                        
                        for tp,t_prime in enumerate(t_primes): 

                            if time_step-t_prime > mem_cutoff:
                                integrand[tp,i,j] = 0
                                continue

                            if t_prime == time_step:
                                integrand[tp,i,j]=eta_t_fit[ts,i,j,0]
                                continue  
                            
                            integrand[tp,i,j] = eta_t_fit[tp,i,j,ts-tp]
                            integrand[tp,i,j] *= velocities_inter[tp,j_atom,j_cart]
                            integrand[tp,j,i] = integrand[tp,i,j] #TODO check integrand[0]
                        force_vec[co,ts,i_atom,i_cart] += np.trapz(integrand[:,i,j],time_axis)

            for ts, time_step in enumerate(inter_time_scale):     
                for i in range(dimension):
                    i_cart = i % 3
                    i_atom = i // 3 
                    
                    nm_work[co,ts] += np.dot(velocities_inter[ts,i_atom,i_cart],force_vec[co,ts,i_atom,i_cart])
                
        self.nm_work = nm_work*self.time_step
        self.force_vec = force_vec

    def calculate_friction_force(self):

        if not hasattr(self,'force_vec'):
            self.frequency_interpolate()
            self.generate_domains()
            self.fourier_transform()
            self.time_interpolate()
            self.convolute()
            self.get_velocities()
            self.velocitiy_interpolation()
            self.mem_integral()
            return(self.force_vec)
        else:
            return(self.force_vec)

    def calculate_work(self):
        import time
        start_time = time.time()
        if not hasattr(self,'nm_work'):
            self.frequency_interpolate()
            print("--- %s FI ---" % (time.time() - start_time))
            self.generate_domains()
            print("--- %s GD ---" % (time.time() - start_time))
            self.fourier_transform()
            print("--- %s FT ---" % (time.time() - start_time))
            self.time_interpolate()
            print("--- %s TI ---" % (time.time() - start_time))
            self.convolute()
            print("--- %s C ---" % (time.time() - start_time))
            self.get_velocities()
            print("--- %s GV ---" % (time.time() - start_time))
            self.velocitiy_interpolation()
            print("--- %s VI ---" % (time.time() - start_time))
            self.mem_integral()
            print("--- %s MI ---" % (time.time() - start_time))
            return(self.nm_work)
        else:
            return(self.nm_work)

        

def Parse_memory_kernels(path_to_calcs,file_range,read=False):

    from coolvib.tools.spectrum import read_memory_kernel
    import numpy as np

    filename  = 'raw_memory.npy'
    bins,re,im,dimension,max_e = read_memory_kernel(path_to_calcs+'/'+str(file_range[0])+'/friction_memory_kernel.out')
    if read:
        print('reading')
        raw_data = np.load(filename)

    else:

        raw_data = np.zeros((len(file_range),dimension,dimension,len(bins)))

        for ts in file_range:        
            path = path_to_calcs+'/'+str(ts)+'/friction_memory_kernel.out'
            try:
                bins,re,im,dimension,max_e = read_memory_kernel(path)
            except:
                print('cannot get mem_kernel for '+str(ts)+' - continuing')
                continue
                    
            raw_data[ts-1,:,:,:] = re
        np.save(filename,raw_data)

    return raw_data,bins,dimension