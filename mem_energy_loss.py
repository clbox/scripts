import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
from ase.db import connect
import math
hbar = _hbar * J * s 
ps = fs*1000

#TODO: Add mode where integrate eta t function to get a 'markov' like tensor and plot energy loss from this
class Postprocessed_memory:

    """ Takes Raw data, array of n_structures,dimension,dimension,energy_bins 
        and calculates the spectral density, Fourier transform, with 
        convolution and calculates energy loss due to friction including
        the effects of memory 
        
        Can do this at multiple cutoff energies if supplied, or one
        INPUT:
        bins / eV
        raw_data / ps-1
        cutoffs / eV
        mem_cutoff / fs
        time_step / fs

        self. all ase units
        
    """
    
    def __init__(self,bins,raw_data,cutoffs,mem_cutoff,friction_indices,time_step,con,debug=False,treat_complex=False):      
        self.raw_data = raw_data / ps
        self.elements= np.shape(raw_data)[1]
        self.cutoffs = cutoffs
        self.mem_cutoff = mem_cutoff * fs
        self.steps = np.shape(raw_data)[0]
        self.friction_indices = friction_indices
        self.time_step = time_step * fs
        self.bins = bins
        self.con = con
        self.dimension = len(friction_indices)*3
        self.debug = debug
        self.treat_complex = treat_complex

    def frequency_interpolate(self):
        """Add extra bins in frequency domain"""
        
        time_bins = np.linspace(4*fs,self.mem_cutoff,500)
        add_bins = (1/(time_bins))*hbar #eV
        new_bins = np.append(self.bins,add_bins)
        new_bins = np.sort(new_bins)

        new_data = np.zeros((self.steps,self.elements,len(new_bins)))

        for ts in range(self.steps):
                for e in range(self.elements):
                    f = interp1d(self.bins,self.raw_data[ts,e,:])
                    new_data[ts,e,:] = f(new_bins)

        self.new_data = new_data
        self.new_bins = new_bins

    def generate_domains(self):
        frequency_list = []
        times_list = []
        self.eta_bar_t_list = []
        for cutoff in self.cutoffs:
            frequencies=((self.new_bins[self.new_bins<cutoff])/hbar) #ase units inverse time
            times = 1/frequencies #ase units time
        
            times = np.append(times,0)
            times = np.sort(times) #check consequencies
            times[-1] = 2**32 #check consequencies assumes 0 freq

            frequency_list.append(frequencies)
            times_list.append(times)
            eta_bar_t = np.zeros((self.steps,self.elements,len(times)))
            self.eta_bar_t_list.append(eta_bar_t)

        self.frequency_list = frequency_list
        self.times_list = times_list

    def get_friction_masses(self):
        """ assumes same atom order for all steps"""
        friction_indices = self.friction_indices
        try:
            atoms = self.con.get_atoms(id=1)
        except:
            print('cannot get atoms for id = '+str(1))
        masses = atoms.get_masses()
        friction_masses = masses[friction_indices]
        return friction_masses

    def fourier_transform(self):
        """Fourier transform to time domain"""
        elements = self.elements
        indices = self.element_index()
        masses = self.get_friction_masses()

        for co in range(len(self.cutoffs)):
            times = self.times_list[co]
            frequencies = self.frequency_list[co]
            func = np.zeros(len(frequencies))
            eta_bar_t = self.eta_bar_t_list[co]
            cos_factor = np.cos(frequencies*times[:,None])

            if self.treat_complex == False:
                for ts in range(self.steps):
                    lambda_omega = self.new_data[ts,:,0:len(frequencies)] #convert from ps-1
                    
                    for e in range(elements):
                        i_atom,j_atom = (indices[e])[0] // 3, (indices[e])[1] // 3
                        func = lambda_omega[e,None,:] * cos_factor * np.sqrt(masses[i_atom]*masses[j_atom])
                        func[:,0]=0
                        eta_bar_t[ts,e,:]=np.trapz(func,frequencies,1)
            # else:
            #TODO: FIGURE THIS OUT
            #     for ts in range(self.steps):
            #         lambda_omega = self.new_data[ts,:,0:len(frequencies)] #convert from ps-1
                    
            #         for e in range(elements):
            #             i_atom,j_atom = (indices[e])[0] // 3, (indices[e])[1] // 3
            #             func = lambda_omega[e,None,:] * cos_factor * np.sqrt(masses[i_atom]*masses[j_atom])
            #             func[:,0]=0
            #             eta_bar_t[ts,e,:]=np.trapz(func,frequencies,1)


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


        dt = (self.time_step/40)
        n_points = math.ceil(self.mem_cutoff/dt)
        time_e_max = dt * (n_points-1)
        n_points = n_points + n_points-1
        times_forward = np.linspace(0,time_e_max,n_points)
        times_up = np.append(times_forward,-times_forward[1:])
        times_up = np.sort(times_up)
        below_zero = np.argwhere(times_up < 0.0)

        for co in range(len(self.cutoffs)):
            
            times = self.times_list[co]   
            eta_bar_t = self.eta_bar_t_list[co]
            
            self.times_up_list.append(times_up)
            
            eta_bar_inter = np.zeros((self.steps,self.elements,len(times_up)))         
            
            for ts in range(self.steps):
                for e in range(self.elements):
                    
                        f = interp1d(times,eta_bar_t[ts,e,:],fill_value="extrapolate")
                        
                        eta_bar_inter[ts,e,:] = f(times_up)
                        
                        eta_bar_inter[ts,e,below_zero] = 0
            

            self.eta_bar_inter_list.append(eta_bar_inter)       


    def convolute(self):
        """Convolute with sinc factor in time domain"""
        self.eta_t_list = []

        for co,cutoff in enumerate(self.cutoffs):
            times_up=self.times_up_list[co]    
            eta_bar_inter = self.eta_bar_inter_list[co]
            dx = times_up[1]-times_up[0]
            #find cutoff

            cutoff_freq = (cutoff/hbar)
            #factors
            sinc = cutoff_freq*np.sinc((times_up*cutoff_freq)/2)    
            
            #exp_factor = np.cos(0.5*times_up*cutoff_freq) + (1j * np.sin(0.5*times_up*cutoff_freq))
            exp_factor = np.cos(0.5*times_up*cutoff_freq)

            #eta_t = np.zeros((final,dimension,dimension,len(times_up)),dtype=complex)
            eta_t = np.zeros((self.steps,self.elements,len(times_up[times_up >= 0.0])))
            
            convolute_factor = sinc*exp_factor
            
            for ts in range(self.steps):
                for e in range(self.elements):
                    eta_t[ts,e,:] = (np.convolve(eta_bar_inter[ts,e,:],convolute_factor,'same')*dx)[times_up >= 0.0]
            self.eta_t_list.append(eta_t)

    def get_velocities(self):
        
        self.all_velocities = np.zeros((self.steps,len(self.friction_indices),3))
        for i in range(self.steps):
            try:
                atoms = self.con.get_atoms(id=i+1)
            except:
                continue
            self.all_velocities[i,:,:] = atoms.get_velocities()[self.friction_indices,:]

    def velocity_interpolation(self):
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
        """Purposefully ensured the interpolated nuclear and 'electronic' time step are equal such that 
        the memory integral can be treated as summation of the reverse diagonals for each nuclear time
        step (then multiplication by the delta t)"""

        old_time_scale = np.linspace(0,(self.steps-1)*self.time_step,self.steps)
        dimension = self.dimension
        inter_time_scale = self.inter_time_scale
        velocities_inter = self.velocities_inter
        force_vec = np.zeros((len(self.cutoffs),len(inter_time_scale),len(self.friction_indices),3))
        nm_work=np.zeros((len(self.cutoffs),len(inter_time_scale)))
        elements = self.elements
        indices = self.element_index()

        dt = inter_time_scale[1]-inter_time_scale[0]

        for co in range(len(self.cutoffs)):
            times_up = self.times_up_list[co]
            eta_t = self.eta_t_list[co]
            times_up = times_up[times_up >= 0.0]
            
            for e in range(elements):
                i = (indices[e])[0]
                j = (indices[e])[1]
                i_cart,j_cart = i % 3, j % 3
                i_atom,j_atom = i // 3, j // 3

                fit = interp1d(old_time_scale,eta_t[:,e,:],kind='linear',axis=0)
                a = fit(inter_time_scale)
                a *= velocities_inter[:,None,j_atom,j_cart] #multiply column wise
                integrand = np.array([np.sum(np.diag(np.fliplr(a), d)) for d in range(np.shape(a)[1]-1,np.shape(a)[1]-np.shape(a)[0]-1, -1)])*dt
                force_vec[co,:,i_atom,i_cart] += integrand
                if i != j:
                    force_vec[co,:,j_atom,j_cart] += integrand

            for i in range(dimension):
                i_cart = i % 3
                i_atom = i // 3
                nm_work[co,:] += velocities_inter[:,i_atom,i_cart]*force_vec[co,:,i_atom,i_cart]

        self.nm_work = nm_work*dt
        self.force_vec = force_vec

    def markov_integral(self):
        """Pseudo markov integral test where eta(t) is integrated for each nuclear time step
        and treated as a markovian friction tensor """

        old_time_scale = np.linspace(0,(self.steps-1)*self.time_step,self.steps)
        dimension = self.dimension
        inter_time_scale = self.inter_time_scale
        velocities_inter = self.velocities_inter
        force_vec = np.zeros((len(self.cutoffs),len(inter_time_scale),len(self.friction_indices),3))
        m_work=np.zeros((len(self.cutoffs),len(inter_time_scale)))
        elements = self.elements
        indices = self.element_index()
        friction_vals = np.zeros((len(self.cutoffs),elements,len(inter_time_scale)))

        dt = inter_time_scale[1]-inter_time_scale[0]

        masses = self.get_friction_masses()

        for co in range(len(self.cutoffs)):
            times_up = self.times_up_list[co]
            eta_t = self.eta_t_list[co]
            times_up = times_up[times_up >= 0.0]
            
            for e in range(elements):
                i = (indices[e])[0]
                j = (indices[e])[1]
                i_cart,j_cart = i % 3, j % 3
                i_atom,j_atom = i // 3, j // 3

                fit = interp1d(old_time_scale,eta_t[:,e,:],kind='linear',axis=0)
                a = fit(inter_time_scale)

                integrand = np.sum(a,axis=1)*dt

                friction_vals[co,e,:] = integrand/np.sqrt(masses[i_atom]*masses[j_atom])

                force_vec[co,:,i_atom,i_cart] += integrand * velocities_inter[:,j_atom,j_cart]
                if i != j:
                    force_vec[co,:,j_atom,j_cart] += integrand

            for i in range(dimension):
                i_cart = i % 3
                i_atom = i // 3
                m_work[co,:] += velocities_inter[:,i_atom,i_cart]*force_vec[co,:,i_atom,i_cart]

        self.m_work = m_work*dt
        self.m_force_vec = force_vec
        self.m_friction_vals = friction_vals

    def calculate_friction_force(self):

        if not hasattr(self,'force_vec'):
            self.frequency_interpolate()
            self.generate_domains()
            self.fourier_transform()
            self.time_interpolate()
            self.convolute()
            self.get_velocities()
            self.velocity_interpolation()
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
            self.velocity_interpolation()
            print("--- %s VI ---" % (time.time() - start_time))
            self.mem_integral()
            print("--- %s MI ---" % (time.time() - start_time))
            return(self.nm_work)
        else:
            return(self.nm_work)

    def calculate_markov_work(self):
        if not hasattr(self,'m_work'):
            self.frequency_interpolate()
            self.generate_domains()
            self.fourier_transform()
            self.time_interpolate()
            self.convolute()
            self.get_velocities()
            self.velocity_interpolation()

            self.markov_integral()
            return(self.m_work)
        else:
            return(self.m_work)

    def element_index(self):
        "Returns list of indices for flattened triangular array"
        a = np.arange(0,self.elements)
        i=0 
        j=0
        indices = []
        for e in a:
            indices.append((i,j))
            if j == self.dimension-1: 
                i += 1 
                j = i 
                continue
            j+=1
        return(indices)
                    

def Parse_memory_kernels(path_to_calcs,file_range,read=False):
    filename  = 'raw_memory.npy'
    bins,re,im,dimension,max_e = read_memory_kernel(path_to_calcs+'/'+str(file_range[0])+'/friction_memory_kernel.out')
    elements = int((((dimension*dimension)-dimension)/2)+dimension)
    if read:
        print('reading')
        raw_data = np.load(filename)

    else:

        raw_data = np.zeros((len(file_range),elements,len(bins)))

        for ts in file_range:        
            path = path_to_calcs+'/'+str(ts)+'/friction_memory_kernel.out'
            try:
                bins,re,im,dimension,max_e = read_memory_kernel(path)
            except:
                print('cannot get mem_kernel for '+str(ts)+' - continuing')
                continue

            if np.max(re) > 10:
                print('Large max for id = ' + str(ts+1) + ' defaulting to zero')
                re=np.zeros((elements,len(bins)))
                im=np.zeros_like(re)
                    
            raw_data[ts-1,:,:] = re
        np.save(filename,raw_data)

    return raw_data,bins

def read_memory_kernel(path,treat_complex=True):
    #In some older calculations the full tensor is printed (old_format = True)
    #In this case need to ignore elements where 
    head_count =0
    old_format = False
    header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction"] #skip lines
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                dimension = int(line.split()[3])
                if int(line.split()[3]) > int(line.split()[4]):
                    continue
                else:
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
    skip = False
    with open(path, "r") as f:
        for line in f:
            if "Friction" in line:
                c=0
                if int(line.split()[3]) > int(line.split()[4]):
                    skip = True
                else:
                    skip = False
                    e +=1
            if any(x in line for x in header):
                continue
            else:
                if skip:
                    continue
                re_memory_kernel[e-1,c]=float(line.split()[1])
                if treat_complex:
                    im_memory_kernel[e-1,c]=float(line.split()[2])
                bins[c]=float(line.split()[0])
                c +=1

    return(bins,re_memory_kernel,im_memory_kernel,dimension,max_e)



class Postprocessed_markov:

    """ Takes connection to database containing each atoms object and 
        friction tensor and calculates the post processed energy loss
        due to (Markovian) electronic friction.
        
        Can do this at multiple cutoff energies if supplied, or one
        INPUT:
        time_step / fs
        database (friction tensor ps^-1)

        Requires the key the friction tensor is stored in the database
        under as a string, e.g 'ft1' 

        Requires the number of steps, the final point in the trajectory
        should be saved in database as id=steps+1

        self. all ase units
        
    """
    
    def __init__(self,steps,friction_indices,time_step,con,key):
        self.friction_indices = friction_indices
        self.time_step = time_step * fs
        self.con = con
        self.dimension = len(friction_indices)*3
        self.key = key
        self.steps = steps
        self.time_scale = np.linspace(0,(self.steps-1)*self.time_step,self.steps)

    def get_all_tensors(self):
        dimension = self.dimension
        steps = self.steps
        all_tensors = np.zeros((steps,dimension,dimension))
        for ts in range(steps):
            row = self.con.get(id=ts+1)
            try:
                tensor = self.string2array(row.get(self.key))
            except:
                'Cannot get tensor for id = ' + str(ts+1) + ' ,defaulting to zero' 
                tensor = np.zeros((dimension,dimension))
            all_tensors[ts,:,:] = tensor / ps
        return all_tensors

    def get_friction_masses(self):
        """ assumes same atom order for all steps"""
        friction_indices = self.friction_indices
        atoms = self.con.get_atoms(id=1)
        masses = atoms.get_masses()
        friction_masses = masses[friction_indices]
        return friction_masses

    def energy_loss(self):
        steps = self.steps
        friction_indices = self.friction_indices
        dimension = self.dimension
        friction_masses = self.get_friction_masses()

        self.m_work = np.zeros((steps))
        self.m_forces = np.zeros((steps,len(friction_indices),3))
        Postprocessed_memory.get_velocities(self)
        all_tensors = self.get_all_tensors() 
    
        for ts in range(steps):
            tensor = all_tensors[ts,:,:]
            vel = self.all_velocities[ts,:,:]
            if np.amax(all_tensors[ts,:,:]) > 8/ps:
                print('id = ' + str(ts) + ' has tensor elements over 8 ps^-1 !, defaulting to zero')
                continue
      
            for i in range(dimension):
                i_atom = i // 3
                for j in range(dimension):
                    j_atom = j // 3
                    mass_factor=np.sqrt(friction_masses[i_atom])*np.sqrt(friction_masses[j_atom])
                    tensor[i,j]=tensor[i,j]*mass_factor 
        
            forces=np.dot(tensor,vel.flatten())
            self.m_work[ts] = np.dot(vel.flatten(),forces)*self.time_step
            self.m_forces[ts,:,:] = forces.reshape(len(friction_indices),3)


    def calculate_friction_force(self):

        if not hasattr(self,'m_forces'):
            self.energy_loss()
            return(self.m_forces)
        else:
            return(self.m_forces)

    def calculate_work(self):

        if not hasattr(self,'m_work'):
            self.energy_loss()
            return(self.m_work)
        else:
            return(self.m_work)

    def string2array(self,string):
        """
        Converts friction tensor that is stored as string in database to
        a numpy array
        """
        dimension = len(string.split(']\n'))
        onedarray = np.fromstring((string.replace('[',' ').replace(']\n',' ')),dtype=float,sep=' ')
        return onedarray.reshape(dimension,dimension)