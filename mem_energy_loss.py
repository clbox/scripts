class Postprocessed_memory:

    """ Takes Raw data, array of n_structures,dimension,dimension,energy_bins 
        and calculates the spectral density, Fourier transform, with 
        convolution and calculates energy loss due to friction including
        the effects of memory 
        
        Can do this at multiple cutoff energies if supplied, or one"""

    def __init__(self,bins,raw_data,cutoffs,mem_cutoff,friction_indices,time_step,con):
        import numpy as np
        from scipy.interpolate import interp1d
        from ase.units import _hbar, J, s, fs
        hbar = _hbar * J
        self.raw_data = raw_data
        self.dimension = np.shape(raw_data)[1]
        self.cutoffs = cutoffs
        self.mem_cutoff = mem_cutoff * fs
        self.steps = np.shape(raw_data)[0]
        self.friction_indices = friction_indices
        self.time_step = time_step * fs
        self.bins = bins

    def frequency_interpolate(self):
        """Add extra bins in frequency domain"""
        
        time_bins = np.linspace(4e-15,self.mem_cutoff/fs,500) #fs
        add_bins = (1/(time_bins))*hbar #eV
        new_bins = np.append(self.bins,add_bins)
        new_bins = np.sort(new_bins)

        new_data = np.zeros((self.steps,self.dimension,self.dimension,len(new_bins)))

        for ts in range(self.steps):
                for i in range(self.dimension):
                    for j in range(self.dimension):
                        if j<i:
                            continue
                        f = interpolate.interp1d(self.bins,self.raw_data[ts,i,j,:])
                        new_data[ts,i,j,:] = f(new_bins)
        self.new_data = new_data
        self.new_bins = new_bins

    def generate_domains(self):
        frequency_list = []
        times_list = []
        self.eta_bar_t_list = []
        for co,cutoff in enumerate(self.cutoffs):
            frequencies=((self.new_bins[self.new_bins<cutoff])/hbar)/s #ase units inverse time
            times = 1/frequencies #ase units time
        
            times = np.append(times,0)

                    
            frequency_list.append(frequencies)
            times_list.append(times)
            eta_bar_t = np.zeros((self.steps,self.dimension,self.dimension,len(times)))
            self.eta_bar_t_list.append(eta_bar_t)

        self.frequency_list = frequency_list
        self.times_list = times_list
    
    def fourier_transform(self,con):
        """Fourier transform to time time domain"""
        from ase import Atoms
        from ase.db import connect

        for ts in range(self.steps):
            try:
                atoms = con.get_atoms(id=ts)
            except:
                print('cannot get atoms for id = '+str(ts+1)+' - continuing')
                continue
            
            masses = atoms.get_masses()
            
            re_memory_kernel=self.new_data[ts,:,:,:]
            
            for i in range(self.dimension):
                i_atom = int(i/3)
                for j in range(self.dimension):
                    if j<i:
                        continue
                    j_atom = int(j/3)
                    mass_factor=np.sqrt(masses[self.friction_indices[i_atom]])*np.sqrt(masses[self.friction_indices[j_atom]]) #amu
                    lambda_omega = re_memory_kernel[i,j,:]*mass_factor/(fs*1000)
            
                    for co in range(len(self.cutoffs)):
                        times = self.times_list[co]
                        frequencies = self.frequency_list[co]
                        func = np.zeros(len(frequencies))
                        eta_bar_t = self.eta_bar_t_list[co]

                        for t in range(len(times)):
                            func=(lambda_omega[0:len(frequencies)])*np.cos(frequencies*times[t])
                            func[0]=0
                            eta_bar_t[ts,i,j,t]=np.trapz(func[:],frequencies[:])
                        self.eta_bar_t_list[co]=eta_bar_t

    def time_interpolate(self):
        """At the moment it is neccesary to interpolate in the  time domain
        as well to get a uniform distribution of bins for the convolution """
        self.times_up_list = []
        self.eta_bar_inter_list = []

        times_up = np.linspace(-self.mem_cutoff*fs,self.mem_cutoff*fs,500)
        below_zero = np.argwhere(times_up < 0.0)
        for co in range(len(self.cutoffs)):
            
            times = self.times_list[co]   
            eta_bar_t = self.eta_bar_t_list[co]
            
            self.times_up_list.append(times_up)
            
            eta_bar_inter = np.zeros((self.steps,self.dimension,self.dimension,len(times_up)))         
            
            for ts in range(self.steps):
                for i in range(self.dimension):       
                    for j in range(self.dimension):
                        if j<i:
                            continue
                    
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
                    for j in range(self.dimension):
                        if j<i:
                            continue  
                        eta_t[ts,i,j,:] = np.convolve(eta_bar_inter[ts,i,j,:],convolute_factor,'same')*dx 
            self.eta_t_list.append(eta_t)

    def get_velocities(self,con):
        
        self.all_velocities = []
        for i in range(self.steps):
            atoms = con.get_atoms(id=i+1)
            self.all_velocities[i,:,:] = atoms.get_velocities()[self.friction_indices,:]

    def velocitiy_interpolation(self):
        """To properly evaluate the memory integral, need to carry out
        integration on a scale denser than the nuclear time step scale"""


        dt = (self.times_up_list[0])[1]-(self.times_up_list[0])[0]
        dts = dt*self.time_step
        #Define time scale as index

        inter_time_scale = np.arange(0,self.steps*time_step,dts)
        old_time_scale = np.arange(0,self.steps*time_step,self.steps)
        self.velocities_inter = np.zeros(inter_time_scale,len(self.friction_indices),3)


        for atom in range(len(self.friction_indices)):
            for cart in range(3):
                linfit = interp1d(old_time_scale,self.all_velocities[:,atom,cart]) 

                self.velocities_inter[:,atom,cart] = linfit(inter_time_scale)

        self.inter_time_scale = inter_time_scale

    def get_atoms(self,con):
        """Just gets the atoms for the first structure to have the masses
            etc to work with"""


    def mem_integral(self,mem_cutoff):
        self.nm_work=np.zeros((len(self.cutoffs),self.steps))
        self.force_vec = np.zeros((len(self.cutoffs),self.steps,len(self.friction_indices),3))
        for co in range(len(self.cutoffs)):

            eta_t = self.eta_t_list[co]*-1
            times = self.times_up_list[co]

            for ts, time_step in enumerate(self.inter_time_scale):

                t_primes = self.inter_time_scale[self.inter_time_scale<=time_step]
                time_axis = time_step - t_primes
                integrand = np.zeros([len(t_primes),self.dimension,self.dimension])

                for tp,t_prime in enumerate(t_primes): 
                    if ts ==0:
                        integrand[0,:,:]=0
                        continue

                    real_below_ts = int((math.floor(t_prime)/2)+1)            
                    real_above_ts = real_below_ts + 1

                    if time_step-t_prime > mem_cutoff:
                        integrand[tp,:,:]=0
                        continue
                    
                    i_cart = -1
                    i_atom = 0
                    for i in range(self.dimension):
                        i_cart = i_cart + 1
                        if (i_cart>2):
                            i_cart = 0
                            i_atom = i_atom + 1
          
                        j_atom = 0
                        j_cart = -1          
                        for j in range(self.dimension):
                            j_cart = j_cart + 1
                            if (j_cart>2):
                                j_cart = 0
                                j_atom = j_atom + 1
                            
                            if j<i:
                                continue

                            eta_t_below = eta_t[real_below_ts-1,i,j,:] 
                            eta_t_above = eta_t[real_above_ts-1,i,j,:] 

                            fit = interp1d([0,1], np.vstack([eta_t_below, eta_t_above]), axis=0) 

                            eta_t_sort = fit(t_prime/((real_above_ts-1)*2))
                            if t_prime == time_step:
                                integrand[tp,i,j]=eta_t_sort[0]
                                continue  

                            integrand[tp,i,j] = inter_bin(time_step-(t_prime*fs),times,eta_t_sort)
                            integrand[tp,i,j]*=self.velocities_inter[tp,j_atom,j_cart]
                            integrand[tp,j,i]= integrand[tp,i,j] #TODO check integrand[0]
                
                i_cart = -1
                i_atom = 0
                for i in range(self.dimension):
                    i_cart = i_cart + 1
                    if (i_cart>2):
                        i_cart = 0
                        i_atom = i_atom + 1
                                
                    j_atom = 0
                    j_cart = -1          
                    for j in range(self.dimension):
                        j_cart = j_cart + 1
                        if (j_cart>2):
                            j_cart = 0
                            j_atom = j_atom + 1

                        self.force_vec[co,ts,i_atom,i_cart]+=np.trapz(integrand[:,i,j],time_axis)

                    self.nm_work[co,ts]+=np.dot(vel[self.friction_indices[i_atom],i_cart],self.force_vec[co,ts,i_atom,i_cart])
                
                
        self.nm_work = self.nm_work*self.time_step

    def calculate_friction_force(self):

        if not hasattr(self,'force_vec'):
            do
        else:
            return(self.force_vec)

    def calculate_work(self):
        
        if not hasattr(self,'nm_work'):
            do
        else:
            return(self.nm_work)

        


    

    
    





    



def Parse_memory_kernels(path_to_calcs,file_range,read=False):

    from coolvib.tools.spectrum import read_memory_kernel

    filename  = 'raw_memory.npy'
    if read:
        print('reading')
        raw_data = np.load(filename)

    else:

        bins,re_memory_kernel,im_memory_kernel,dimension,max_e = read_memory_kernel(path_to_calcs+str(file_range[0])+'/friction_memory_kernel.out')

        raw_data = np.zeros((file_range[-1]-1,dimension,dimension,len(bins)))
        for ts in range(file_range):
            row = con.get(id=ts)
            try:
                atoms = con.get_atoms(id=ts)
            except:
                print('cannot get atoms for '+str(ts)+' - continuing')
                continue
        
            path = path_to_calcs+str(ts)+'/friction_memory_kernel.out'
            try:
                bins,re_memory_kernel,im_memory_kernel,dimension,max_e = read_memory_kernel(path)
            except:
                print('cannot get mem_kernel for '+str(ts)+' - continuing')
                continue
                    
            raw_data[ts-1,:,:,:] = re_memory_kernel
        np.save(filename,raw_data)

        return raw_data