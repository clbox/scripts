import numpy as np
import sys
import os
from ase import Atoms
from scipy import special
from ase.units import _hbar, J, s, fs
from scipy.integrate import simps
import time
start_time = time.time()
hbar = _hbar * J * s
ps = fs*1000
boltzmann_kB = 8.617333262145e-5 #eV K-1
hartree = 27.2113845
#TODO spin

def fermi_pop(x,x0,T):
    if T<1e-10:
        fermi_pop = np.zeros_like(x)
        for i,ex in enumerate(x):
            if ex<x0:
                fermi_pop[i] = 1.0
            else:
                fermi_pop[i] = 0.0
    else:
        fermi_pop = (1./(np.exp((x-x0)/(boltzmann_kB*T))+1))
    return fermi_pop

def gaussian_function(x,x0,s):
    return 1./np.sqrt(2*np.pi*s*s)*np.exp(-0.5*(((x-x0) / (s))**2))

def gaussian_function2(x,x0,s):
    return 1./np.sqrt(np.pi*s*s)*np.exp(-(((x-x0) / (s))**2))

def gaussian_function_qe(ei,ej,ef,s):

    return 1./(np.pi*s**2) * np.exp(-((ef-np.array(ei))**2 + (ef-np.array(ej))**2)/s**2)

# Franceso mauri 2005
def gaussian_function_mauri(x,x0,s): 
    return 1./np.sqrt(np.pi*s)*np.exp((((x) / (s))**2))

def gaussian_norm(x0, s):
    return 0.5 * (1-special.erf((-x0/s)*(1/np.sqrt(2))))

def gaussian_norm2(x0, s):
    return 0.5 * (1-special.erf((-x0)))

def gaussian_derivative(x,x0,s):

    gd  = (-x+x0)/(np.sqrt(2*np.pi)*s**3) * np.exp((-(x-x0)**2)/(2*s**2))

    return gd

def fermi_derivative(x,x0,T):

    fermi_derivative = np.exp((x-x0)/(boltzmann_kB*T))
    fermi_derivative = fermi_derivative / (((fermi_derivative+1)**2)*T*boltzmann_kB)
        
    return fermi_derivative

def lorentzian_function(x,x0,s):
    return (1./np.pi)*((0.5*s)/((x-x0)*(x-x0)+(0.5*s)*(0.5*s)))

def methfessel_paxton_function(e,e0,s,N):
    x = (e-e0)/s 

    delta_function_methfessel = 0.

    for i in range(0, N+1):
        tmp = 0.
        if i==0:
            tmp = np.exp(-x*x)/np.sqrt(np.pi)
        elif i==1:
            tmp = -(4.*x*x-2)*np.exp(-x*x)/(4*np.sqrt(np.pi))
        elif i==2:
            tmp = (16.*x**4-48.*x*x*12)*np.exp(-x*x)/(32.*np.sqrt(np.pi))
        elif i==3:
            tmp = -(64.*x**6-480.*x**4+720.*x*x-120.)*np.exp(-x*x)/(384.*np.sqrt(np.pi))
        elif i==4:
            tmp = (256.*x**8-3584.*x**6+13440.*x**4-13440.*x*x+1680.)*\
                np.exp(-x*x)/(6144.*np.sqrt(np.pi))
        elif i==5:
            tmp = -(1024.*x**10-23040.*x**8+161280.*x**6-403200.*x**4+302400.*x*x-30240.)*\
                np.exp(-x*x)/(122880.*np.sqrt(np.pi))
        elif i==6:
            tmp = (1024.*x**10-23040.*x**8+161280.*x**6-403200.*x**4+302400.*x*x-30240.)*\
                np.exp(-x*x)/(122880.*np.sqrt(np.pi))
        else:
            print('delta_function_methfessel: N>6 not allowed')
        delta_function_methfessel = delta_function_methfessel + tmp
    return delta_function_methfessel

# Changed 2023
class friction_gamma_parser():
    def __init__(self,aims_file,gamma_files,nspin):
        self.nspin = nspin
        self.chem_pot = self.parse_chem_pot(aims_file,self.nspin)

        self.chem_pot_response = self.parse_chem_pot_response(aims_file,)
        #gamma_files.sort()
        max_k,gamma_files_sorted = self.sort_gamma_files(gamma_files)
        self.friction_masses = self.parse_fiction_masses(aims_file)
        eis,ejs,couplings,kweights = self.parse_gamma_couplings(gamma_files_sorted)




        self.ndim = len(self.friction_masses) * 3


        #self.ks = np.array(ks)
        self.max_k = max_k
        #self.coords = np.array(coords)
        self.eis = eis
        self.ejs = ejs
        self.couplings = couplings
        self.kweights = np.array(kweights)



    def sort_gamma_files(self,gamma_files):

        #gamma_files_sorted = np.array(gamma_files)

        
        if any('spin_002' in item for item in gamma_files):
            self.nspin = 2
        else:
            self.nspin =1


        if self.nspin == 2:
            assert(len(gamma_files)%2==0)
        
        gamma_files_sorted = np.empty(shape=(int(len(gamma_files)/self.nspin), self.nspin), dtype=object)


        if self.nspin ==2:
            n1=0
            n2=0
            for i, gf in enumerate(gamma_files): 
                if "_spin_001" in gf:
                    gamma_files_sorted[n1,0] = gf
                    n1+=1
                elif "_spin_002" in gf:
                    gamma_files_sorted[n2,1] = gf
                    n2+=1
        else:
            gamma_files_sorted[:,0] = gamma_files

        k_no_list = []
        for i, gf in enumerate(gamma_files_sorted[:,0]):
            if "spin" in gf:
                k_no = int(gf.replace('friction_gamma_k','').replace('.out','').replace("_spin_001",""))
            else:
                k_no = int(gf.replace('friction_gamma_k','').replace('.out',''))

            k_no_list.append(k_no)

        idx = np.argsort(k_no_list)

        gamma_files_sorted[:,:] = gamma_files_sorted[idx,:]
        k_no_list = np.array(k_no_list)[idx]

        max_k = np.max(k_no_list)


        return max_k,gamma_files_sorted

    def parse_chem_pot(self,aims_file,nspin):
        chem_pot = np.zeros((2))
        with open(aims_file, "r") as af:
            for line in af:
                if '**FRICTION**' in line:
                    break
                if '| Chemical potential (Fermi level):' in line:
                    chem_pot[:] = float(line.split()[-2])
                if '| Chemical potential, spin up:' in line:
                    chem_pot[0] = float(line.split()[-2]) 
                if '| Chemical potential, spin dn:' in line:
                    chem_pot[1] = float(line.split()[-2])  
        return chem_pot

    def parse_chem_pot_response(self,aims_file):
        chem_pot_response = []
        with open(aims_file, "r") as af:
            for line in af:
                if ' | Response of chemical potential (Fermi level):' in line:
                    chem_pot_response.append(float(line.split()[-2]))

        chem_pot_response = np.array(chem_pot_response)
        
        return chem_pot_response

    def parse_fiction_masses(self,aims_file):
        friction_atoms = []
        with open(aims_file, "r") as af:
            read = False
            for line in af:
                if 'The contents of geometry.in will be repeated verbatim below' in line:
                    read=True
                if "moment" in line:
                    continue
                if 'calculate_friction .true.' in line:
                    if (len(element)>2):
                        continue
                    else:
                        friction_atoms.append(element)
                if read:
                    try:
                        element = line.split()[-1]
                    except:
                        continue
        
        a = Atoms(symbols=friction_atoms)
        friction_masses = a.get_masses()

        return friction_masses

    def parse_gamma_couplings(self,gamma_files_sorted):
        #Parsing a limited amount to reduce memory usage
        print("--- %s Start parser ---" % (time.time() - start_time))

        # eis = [[0 for columns in range(0)] for i_rows in range(len(gamma_files))]

        #ejs = [[0 for columns in range(0)] for i_rows in range(len(gamma_files))]



        couplings = [[0 for columns in range(0)] for i_rows in range(self.nspin)]
        eis = [[0 for columns in range(0)] for i_rows in range(self.nspin)] 
        ejs = [[0 for columns in range(0)] for i_rows in range(self.nspin)]  


        for i_spin in range(self.nspin):
            couplings_s = [[0 for columns in range(0)] for i_rows in range(len(gamma_files_sorted[:,0]))]
            eis_s = [[] for i_rows in range(len(gamma_files_sorted[:,0]))]
            ejs_s = [[] for i_rows in range(len(gamma_files_sorted[:,0]))]

            kweights = []
            ks = [] 
            for i,file in enumerate(gamma_files_sorted[:,i_spin]):
                #print(file)
                couplings_k = []
                eis_k = []
                ejs_k = []

                with open(file, "r") as f:
                    read_gamma = False

                    for line in f:
                        if '| k-point:' in line:
                            #k = int(line.split()[2])
                            kweight = float(line.split()[-1])
                            kweights.append(kweight)
                        elif 'Friction component ' in line:
                            read_gamma = True
                            coord = int(line.split()[-1])-1
                            #coords.append(coord)
                        elif read_gamma:
                            
                            ei = float(line.split()[0])
                            ej = float(line.split()[1])
                            
                            # Trim eigenvalues
                            if ei>ej:
                                continue

                            if ei == ej:
                            #     # if abs(ei-self.chem_pot) > 0.05:
                            #     if not (ei==self.chem_pot):
                               continue


                            coupling = complex(float(line.split()[2]),float(line.split()[3])) 
                            couplings_k.append(coupling)

                            if coord==1: #Eigenvalues same for each coomponent
                                eis_k.append(ei)
                                ejs_k.append(ej)

                    couplings_s[i] = couplings_k
                    eis_s[i] = eis_k
                    ejs_s[i] = ejs_k

            couplings[i_spin] = couplings_s
            eis[i_spin] = eis_s
            ejs[i_spin] = ejs_s

                        
        print("--- %s End parser ---" % (time.time() - start_time))
        return eis,ejs,couplings,kweights

class friction_tensor():
    #So ks are base zero (i.e first k_point = 0)
    #ks is just a list of what k point each excitation is associated with
    def __init__(self,parser,temp,sigma,nspin):
        self.couplings = parser.couplings
        # self.coords = parser.coords
        self.eis = parser.eis
        self.ejs = parser.ejs
        #self.ks = parser.ks
        self.max_k = parser.max_k
        self.chem_pot = parser.chem_pot
        self.temp = temp
        self.kweights=parser.kweights
        self.sigma = sigma
        self.friction_masses = parser.friction_masses
        self.nspin = nspin

        self.ndim = parser.ndim

    def calc_tensor(self,mode='default',order=2,perturbing_energy=0.,fermi_correction='none',chemical_potential=None):
        print("--- %s Start calc tensor ---" % (time.time() - start_time))
        #ks = self.ks
        # coords = self.coords
        # ndim = np.max(coords)+1

        ndim = self.ndim
        #max_k = np.max(ks)
        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings
        tensor = np.zeros((ndim,ndim))

        if chemical_potential is None:
            chem_pot = self.chem_pot
        else:
            chem_pot = chemical_potential


        temp = self.temp
        nspin = self.nspin

        if perturbing_energy > 0:
            if mode in ['default','double_delta','prb_print','double_delta_methfessel_paxton']:
                raise ValueError('mode can only deal with 0 perturbing energy')
            

        for i_spin in range(self.nspin):

            couplings_s = couplings[i_spin]
            #Loop over k_points
            
            for i_k,k in enumerate(range(1,self.max_k+1)):

                kw = self.kweights[i_k]

                es = np.array(ejs[i_spin][i_k])-np.array(eis[i_spin][i_k])


                if mode == 'no_norm':
                    gaussian_delta_factor = ((gaussian_function(es,perturbing_energy,self.sigma))/(perturbing_energy))*kw*2/nspin 
                else:
                    gaussian_delta_factor = ((gaussian_function(es,perturbing_energy,self.sigma))/(es))*kw*2/nspin


                fermi_factor = fermi_pop(np.array(eis[i_spin][i_k]),chem_pot[i_spin],temp)-fermi_pop(np.array(ejs[i_spin][i_k]),chem_pot[i_spin],temp)


                couplings_k = np.array_split(np.array(couplings_s[i_k]), ndim)

                omegas = es/hbar

                for i in range(ndim):
                    #i_idx = np.where((coords == i) & (ks == k))[0]
                    for j in range(ndim):
                        if j < i:
                            continue
                        #j_idx = np.where((coords == j) & (ks == k))[0]

                        nacs = np.conjugate(couplings_k[i])*couplings_k[j]

                        # for c_t,coupling in enumerate(nacs):
                        #     if eis[i_k][c_t] == eis[i_k][c_t]:
                        #         # nacs[c_t] = np.real(nacs[c_t]) + 0j
                        #         nacs[c_t] = np.conjugate((couplings_k[i][c_t] - e_f_changes[i])) * (couplings_k[j][c_t] - e_f_changes[j])


                        if mode=='default': # With additional normalisation described in RJM PRB 2016
                            # We no longer like this
                    
                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor*gaussian_delta_factor/gaussian_norm(es,self.sigma))

                        elif mode=='double_delta':
                
                            tensor[i,j] += np.sum(nacs*\
                                (gaussian_function(np.array(eis[i_spin][i_k])-chem_pot[i_spin],0,self.sigma)*gaussian_function(np.array(ejs[i_spin][i_k])-chem_pot[i_spin],0,self.sigma))*\
                                kw*2/nspin)




                        elif mode=='occ_delta_only':
                
                            tensor[i,j] += np.sum(
                                (gaussian_function(np.array(ejs[i_spin][i_k])-chem_pot[i_spin],0,self.sigma))*\
                                kw*2/nspin)


                        elif mode=='unocc_delta_only':
                
                            tensor[i,j] += np.sum(gaussian_function(np.array(ejs[i_spin][i_k])-chem_pot[i_spin],0,self.sigma)*\
                                kw*2/nspin)


                        elif mode=='double_delta_mauri':
                
                            tensor[i,j] += np.sum(nacs*\
                                (gaussian_function_mauri(np.array(eis[i_spin][i_k]),0,self.sigma)*gaussian_function_mauri(np.array(ejs[i_spin][i_k]),0,self.sigma))*\
                                kw*2/nspin)


                        elif mode=='novko':
                
                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor*(gaussian_derivative(np.array(eis[i_spin][i_k]),np.array(ejs[i_spin][i_k]),self.sigma))*\
                                kw*2/nspin)

                        elif mode=='double_delta_qe':

                                tensor[i,j] += np.sum(nacs*\
                                (gaussian_function_qe(eis[i_spin][i_k],ejs[i_spin][i_k],chem_pot[i_spin],self.sigma))* \
                                kw*2/nspin)

                        elif mode=='prb_print':
                
                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor/(es)\
                                *(gaussian_function(es,0,self.sigma)/gaussian_norm2(es,self.sigma))*kw*2/nspin)


                        elif mode=='no_norm':

                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor*gaussian_delta_factor)

                        elif mode=='no_norm_approx':

                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor*gaussian_delta_factor)


                        elif mode=='fermi_derivative':
            
                            tensor[i,j] += np.sum(nacs*\
                                fermi_derivative(es,chem_pot[i_spin],temp)\
                                *(gaussian_function(es,0,self.sigma))*kw*2/nspin)

                        elif mode == 'laplace':

                            omega_factors = (((perturbing_energy/hbar)*omegas)/(((perturbing_energy/hbar)**2)+omegas**2))
                            tensor[i,j] += np.sum((nacs/(es**2))*fermi_factor*omega_factors*kw*2/nspin)

                        elif mode=='gaussian2_no_norm':

                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor/(es)\
                                *(gaussian_function2(es,perturbing_energy,self.sigma))*kw*2/nspin)
                                
                        elif mode=='lorentzian': #no additional normalisation
                
                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor/(es)\
                                *(lorentzian_function(es,perturbing_energy,self.sigma))*kw*2/nspin)

                        elif mode=='methfessel_paxton': #no additional normalisation
                    
                            tensor[i,j] += np.sum(nacs*\
                                fermi_factor/(es)\
                                *(methfessel_paxton_function(es,perturbing_energy,self.sigma,order))*kw*2/nspin)

                        elif mode=='double_delta_methfessel_paxton':
                    
                            tensor[i,j] += np.sum(nacs*\
                                (methfessel_paxton_function(np.array(ejs[i_spin][i_k])-chem_pot[i_spin],0,self.sigma,order)*methfessel_paxton_function(np.array(ejs[i_spin][i_k])-chem_pot[i_spin],0,self.sigma,order))*\
                                kw*2/nspin)

                        else:
                            print('No viable tensor mode selected')
                        
        #tensor *= hbar*np.pi*2 
        #CLB 2021: don't think should have factor 2, removing it gives right tensor value
        tensor *= hbar*np.pi

        for i in range(ndim):
            mass_i = self.friction_masses[i // 3]
            for j in range(ndim):
                if j < i:
                    continue
                mass_j = self.friction_masses[j // 3]
                tensor[i,j] = tensor[i,j]/np.sqrt(mass_i*mass_j)
                tensor[j,i] = tensor[i,j]

        print("--- %s End calc tensor ---" % (time.time() - start_time))
        return tensor*ps


class laplacian():
    #So ks are base zero (i.e first k_point = 0)
    #ks is just a list of what k point each excitation is associated with
    def __init__(self,parser,temp,nspin):
        self.couplings = parser.couplings
        self.coords = parser.coords
        self.eis = parser.eis
        self.ejs = parser.ejs
        #self.ks = parser.ks
        self.max_k = parser.max_k
        self.chem_pot = parser.chem_pot
        self.temp = temp
        self.kweights=parser.kweights
        self.friction_masses = parser.friction_masses
        self.nspin = nspin

    def calc_laplacian(self,mode='default',sigma=0.01,min_e=0.01,max_e=1.0,chemical_potential=None):
        print("--- %s Start calc tensor ---" % (time.time() - start_time))
        #ks = self.ks
        coords = self.coords
        ndim = np.max(coords)+1
        #max_k = np.max(ks)
        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings


        if chemical_potential is None:
            chem_pot = self.chem_pot
        else:
            chem_pot = chemical_potential

        temp = self.temp
        nspin = self.nspin

        e_grid = np.linspace(min_e,max_e,100)
        discretization_length = e_grid[1]-e_grid[0]

        spectrum = np.zeros(len(e_grid))

        laplacian = np.zeros((ndim,ndim,len(e_grid)))

        #Loop over k_points 
        for i_k,k in enumerate(range(1,self.max_k+1)):

            kw = self.kweights[k-1]
            es = np.array(ejs[i_k])-np.array(eis[i_k])


            fermi_factor = (fermi_pop(np.array(eis[i_k]),chem_pot,temp)-fermi_pop(np.array(ejs[i_k]),chem_pot,temp))*kw*2/nspin 

            couplings_k = np.array_split(np.array(couplings[i_k]), ndim)

            omegas = es/hbar



            for i in range(ndim):
                #i_idx = np.where((coords == i) & (ks == k))[0]
                for j in range(ndim):
                    if j < i:
                        continue
                    #j_idx = np.where((coords == j) & (ks == k))[0]

                    nacs = np.conjugate(couplings_k[i])*couplings_k[j]


                    # for c_t,coupling in enumerate(nacs):
                    #     if eis[i_k][c_t] == eis[i_k][c_t]:
                    #         # nacs[c_t] = np.real(nacs[c_t]) + 0j
                    #         nacs[c_t] = np.conjugate((couplings_k[i][c_t] - e_f_changes[i])) * (couplings_k[j][c_t] - e_f_changes[j])


                    if mode=='default': 

                        for e_counter,e in enumerate(e_grid):
                            omega_factors = (((e_grid[e_counter]/hbar)*omegas)/(((e_grid[e_counter]/hbar)**2)+omegas**2))
                            laplacian[i,j,e_counter] += np.sum((nacs/(es**2))*fermi_factor*omega_factors)


                    elif mode == 'gaussian': #TODO: NUMERICAL NORMALISATION OF GAUSSIAN?

                        # for c_t,coupling in enumerate(nacs):

                            gaussian_norm = 0.
                            spectrum[:] = 0.

                            gaussian_contribution = np.zeros((len(es)))

                            for e_counter,e in enumerate(e_grid):

                                gaussian_contribution[:] += gaussian_function(e,es,sigma)


                            for e_counter,e in enumerate(e_grid):

                                # gaussian_norm += gaussian_function(e,es,sigma)
                                # laplacian[i,j,e_counter] += np.sum((nacs[c_t])*(fermi_factor[c_t]/es[c_t])*gaussian_function(es[c_t],e,sigma))
                                spectrum[e_counter] = np.sum((nacs)*(fermi_factor/es)*gaussian_function(e,es,sigma)/(gaussian_contribution[:]*discretization_length))


                            laplacian[i,j,:] += spectrum
            

                    else:
                        print('No viable tensor mode selected')
                    
        #tensor *= hbar*np.pi*2 
        #CLB 2021: don't think should have factor 2, removing it gives right tensor value
        laplacian *= hbar*np.pi

        for i in range(ndim):
            mass_i = self.friction_masses[i // 3]
            for j in range(ndim):
                if j < i:
                    continue
                mass_j = self.friction_masses[j // 3]
                laplacian[i,j,:] = laplacian[i,j,:]/np.sqrt(mass_i*mass_j)
                laplacian[j,i,:] = laplacian[i,j,:]

        print("--- %s End calc tensor ---" % (time.time() - start_time))
        return e_grid,laplacian*ps

class calc_gamma():
    #discretize ei and ejs with gaussian functions of smearing width sigma
    def __init__(self,ks,coords,eis,ejs,couplings,kweights,chem_pot,sigma,nspin,min_e,max_e,npoints):
        self.couplings = couplings
        self.coords = coords
        self.eis = eis
        self.ejs = ejs
        self.ks = ks
        self.chem_pot = chem_pot
        self.kweights=kweights
        self.sigma = sigma
        self.nspin = nspin
        self.min_e = min_e
        self.max_e = max_e
        self.npoints = npoints

    def dos_binning(self):
        print("--- %s Start dos binning ---" % (time.time() - start_time))
        min_e = self.min_e
        max_e = self.max_e
        npoints = self.npoints
        coords = self.coords
        ndim = np.max(coords)

        h_grid = np.linspace(min_e,max_e,npoints) #eV
        e_grid =  np.linspace(min_e,max_e,npoints)

        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings


        gamma = np.zeros((ndim+1,len(e_grid),len(h_grid)),dtype=np.complex128)

        for n in range(ndim+1):
            idx = np.where((coords == n))[0]
            cs = couplings[idx]
            eis_n = eis[idx]
            ejs_n = ejs[idx]
            for i,c in enumerate(cs):
                gauss_mat = np.outer(gaussian_function(eis_n[i],e_grid,0.01),gaussian_function(ejs_n[i],h_grid,0.01))
                gamma[n,:,:] +=gauss_mat*c

        print("--- %s End dos binning ---" % (time.time() - start_time))
        return gamma

#class calc_time_tensor():
    # #discretize ei and ejs with gaussian functions of smearing width sigma
    # def __init__(self,ks,coords,eis,ejs,couplings,kweights,chem_pot,nspin,min_e,max_e,npoints,friction_masses,temp=300,sigma=0.01):
    #     self.couplings = couplings
    #     self.coords = coords
    #     self.eis = eis
    #     self.ejs = ejs
    #     self.ks = ks
    #     self.chem_pot = chem_pot
    #     self.kweights=kweights
    #     self.temp = temp
    #     self.sigma = sigma
    #     self.nspin = nspin
    #     self.min_e = min_e
    #     self.max_e = max_e
    #     self.npoints = npoints
    #     self.ndim = np.max(coords)
    #     self.friction_masses = friction_masses

    # def get_e_grid(self):
    #     min_e = self.min_e
    #     max_e = self.max_e
    #     npoints = self.npoints
    #     e_grid =  np.linspace(min_e,max_e,npoints)
    #     return e_grid

    # def get_h_grid(self):
    #     min_e = self.min_e
    #     max_e = self.max_e
    #     npoints = self.npoints
    #     h_grid = np.linspace(min_e,max_e,npoints)
    #     return h_grid
    
    # def get_ex_grid(self):
    #     min_e = self.min_e
    #     max_e = self.max_e
    #     npoints = self.npoints
    #     h_grid = np.linspace(min_e,max_e,npoints)
    #     return h_grid

    # def calc_A(self):
    #     print("--- %s Start calc A ---" % (time.time() - start_time))

    #     coords = self.coords
    #     ndim = self.ndim
    #     temp = self.temp
    #     chem_pot = self.chem_pot
    #     e_grid = self.get_e_grid()
    #     h_grid = self.get_h_grid()
    #     eis = self.eis
    #     ejs = self.ejs
    #     couplings = self.couplings
    #     sigma = self.sigma

    #     A = np.zeros((ndim+1,ndim+1,len(e_grid),len(h_grid)),dtype=np.complex128)

    #     for n1 in range(ndim+1):
    #         for n2 in range(ndim+1):
    #             if n2< n1:
    #                 continue
    #             idx = np.where((coords == n1))[0]
    #             idx2 = np.where((coords == n2))[0]
    #             cs = couplings[idx]
    #             cs2 = couplings[idx2]
    #             eis_n = eis[idx]
    #             ejs_n = ejs[idx2]
    #             for i,c2 in enumerate(cs2):
    #                 c = np.conjugate(cs[i])
    #                 gauss_mat = np.outer(gaussian_function(eis_n[i],e_grid,sigma),gaussian_function(ejs_n[i],h_grid,sigma))\
    #                     *(fermi_pop(eis_n[i],chem_pot,temp)-fermi_pop(ejs_n[i],chem_pot,temp))
    #                 #print(np.max(gaussian_function(ejs_n[i],h_grid,0.01)))
    #                 A[n1,n2,:,:] +=gauss_mat*c*c2
    #     print("--- %s End calc A ---" % (time.time() - start_time))
    #     return A
    
    # def calc_A2(self):
    #     print("--- %s Start calc A2 ---" % (time.time() - start_time))
    #     coords = self.coords
    #     ndim = self.ndim
    #     temp = self.temp
    #     chem_pot = self.chem_pot
    #     ex_grid = self.get_ex_grid()
    #     eis = self.eis
    #     ejs = self.ejs
    #     couplings = self.couplings
    #     sigma = self.sigma
    #     de = ex_grid[1]-ex_grid[0]

    #     A = np.zeros((ndim+1,ndim+1,len(ex_grid)),dtype=np.complex128)

    #     for n1 in range(ndim+1):
    #         idx = np.where((coords == n1))[0]
    #         cs = couplings[idx]
    #         eis_n = eis[idx]
    #         for n2 in range(ndim+1):
    #             if n2< n1:
    #                 continue
    #             idx2 = np.where((coords == n2))[0]
    #             cs2 = couplings[idx2]
    #             ejs_n = ejs[idx2]
    #             for i,c2 in enumerate(cs2):
    #                 c = np.conjugate(cs[i])*c2
    #                 gauss_mat = gaussian_function(ejs_n[i]-eis_n[i],ex_grid,sigma)*gaussian_function(ejs_n[i]-eis_n[i],ex_grid,sigma)
    #                     #*(fermi_pop(eis_n[i],chem_pot,temp)-fermi_pop(ejs_n[i],chem_pot,temp))
    #                 A[n1,n2,:] +=gauss_mat*c/(np.sum(gauss_mat)*de)
    #     print("--- %s End calc A2 ---" % (time.time() - start_time))
    #     return A

    # def mass_weight(self,tensor):
    #     ndim = self.ndim
    #     for i in range(ndim+1):
    #         mass_i = self.friction_masses[i // 3]
    #         for j in range(ndim+1):
    #             if j < i:
    #                 continue
    #             mass_j = self.friction_masses[j // 3]
    #             tensor[i,j] = tensor[i,j]/np.sqrt(mass_i*mass_j)
    #             tensor[j,i] = tensor[i,j]
    #     return tensor

    # def evaluate_tensor(self):
    #     print("--- %s Start eval tensor ---" % (time.time() - start_time))
    #     ndim = self.ndim
    #     h_grid = self.get_h_grid()
    #     e_grid =  self.get_e_grid()
    #     A = self.calc_A()
    #     for i in range(len(e_grid)):
    #         for j in range(len(h_grid)):
    #             if (h_grid[j]-e_grid[i])==0:
    #                 A[:,:,i,j]=0
    #             else:
    #                 A[:,:,i,j] = A[:,:,i,j]/(h_grid[j]-e_grid[i])


    #     tensor = np.zeros((ndim+1,ndim+1))
    #     for n1 in range(ndim+1):
    #         for n2 in range(ndim+1):
    #             if n2< n1:
    #                 continue
    #             b = simps(A[n1,n2,:,:],e_grid,axis=0)
    #             tensor[n1,n2] = simps(b,h_grid)
    #             tensor[n2,n1] = tensor[n1,n2]

    #     tensor = self.mass_weight(tensor)
    #     tensor *= 2*ps*ps
    #     print("--- %s End eval tensor ---" % (time.time() - start_time))
        return tensor

    
    def evaluate_tensor2(self):
        print("--- %s Start eval tensor2 ---" % (time.time() - start_time))
        ndim = self.ndim
        ex_grid = self.get_ex_grid()
        A = self.calc_A2()
        for i in range(len(ex_grid)):
            if (ex_grid[i])==0:
                A[:,:,i]=0
            else:
                A[:,:,i] = A[:,:,i]/(ex_grid[i])


        tensor = np.zeros((ndim+1,ndim+1))
        for n1 in range(ndim+1):
            for n2 in range(ndim+1):
                if n2< n1:
                    continue
                tensor[n1,n2] = simps(A[n1,n2,:],ex_grid)
                tensor[n2,n1] = tensor[n1,n2]

        tensor = self.mass_weight(tensor)
        tensor *= 2*ps*ps
        print("--- %s End eval tensor2 ---" % (time.time() - start_time))
        return tensor

# ..... 2021 rewrite 
class friction_output_parser_2021():
    def __init__(self):
        print('Initialized')

    def parse_chem_pot(self,aims_file):
        chem_pot = 0
        with open(aims_file, "r") as af:
            for line in af:
                if '--FRICTION--' in line or '**FRICTION**' in line:
                    break
                if '| Chemical potential (Fermi level):' in line:
                    chem_pot = float(line.split()[-2])
        return chem_pot




    def parse_fiction_masses(self,aims_file):
        friction_atoms = []
        with open(aims_file, "r") as af:
            read = False
            for line in af:
                if 'The contents of geometry.in will be repeated verbatim below' in line:
                    read=True
                if 'calculate_friction .true.' in line:
                    friction_atoms.append(element)
                if read:
                    try:
                        element = line.split()[-1]
                    except:
                        continue
        
        a = Atoms(symbols=friction_atoms)
        friction_masses = a.get_masses()

        return friction_masses

    def parse_gamma_couplings(self,gamma_files):
        #Parsing all for now but in future less memory intensive to 
        #process whilst parsing
        print("--- %s Start parser ---" % (time.time() - start_time))

        print('Sort gamma files')
        gamma_files.sort()


        k_points = np.zeros((len(gamma_files),4)) # x  y z weight
        couplings = []
        eigenvalues = []
        n_excits = []


        for i,file in enumerate(gamma_files):
            #print(file)
            with open(file, "r") as f:
                read_gamma = False
                for line in f:
                    if '| k-point:' in line:
                        k_points[i,0] = float(line.split()[4]) # x
                        k_points[i,1] = float(line.split()[5]) # y
                        k_points[i,2] = float(line.split()[6]) # z
                        k_points[i,3] = float(line.split()[-1]) # k_weight
                        continue

                    if 'Friction component for' in line:
                        dimension = int(line.split()[-1])
                        continue


            #  n_excits, 4
            # 0 = ei , 1 = ej , 2 = Re(Y) , 3 = Im(Y)
            gamma_k_data = np.loadtxt(file,comments='Friction',skiprows=4)
            n_excits_k = int(np.shape(gamma_k_data)[0] / dimension)

            gamma_k_reshape = np.reshape(gamma_k_data,(n_excits_k,dimension,4))

            eigenvalues_k = np.zeros((n_excits_k,2))
            eigenvalues_k[:,:] = gamma_k_reshape[:,0,:2] #only  need first dimension since idetical for each coordinate [ground state]

            couplings_k = np.zeros((n_excits_k,dimension),dtype=np.cdouble)
            couplings_k = gamma_k_reshape[:,:,2] + 1j * gamma_k_reshape[:,:,3] 
            

            couplings.append(couplings_k)
            eigenvalues.append(eigenvalues_k)  
            n_excits.append(n_excits_k)
            
        print("--- %s End parser ---" % (time.time() - start_time))

        return eigenvalues, couplings, n_excits, k_points
                
class calc_time_tensor_2021():

    def __init__(self,chem_pot,e_cutoff,friction_masses,temp=300,sigma=0.01):

        # couplings  [2,n_spin,n_k_points,dimension,3,max_n_excits]
        # coupling list (n_k_points long) of arrays of dimension n_excits,dimension

        # first dimension is for each structure
        # penultimate dimension:  0 = ei,  1 = ej , 2 = coupling


        # self.eigenvalues = eigenvalues
        # self.couplings = couplings
        # self.n_excits = n_excits


        #self.n_spin = np.shape(couplings)[1] #number of spin dimensions, 1 or 2
        self.n_spin = 1


        # self.n_k_points = np.shape(couplings)[2]
        # self.n_k_points = len(couplings)


        #self.n_dim = np.shape(couplings)[3] # dimension of tensor = friction_atoms * 3
        #self.n_dim = np.shape(couplings[0])[1]
        
        #self.max_n_excits = np.shape(couplings)[-1]

        self.chem_pot = chem_pot # Chemical potential (i.e Fermi level) / eV
        #self.k_weights=k_weights # list of k_weights , same order as n_k_points dimension
        self.temp = temp  # Temperature / K
        self.sigma = sigma # Broadening / eV

        self.e_cutoff = e_cutoff # energy cutoff / eV

        self.friction_masses = friction_masses # Masses of friction atoms, in same order as dimension (/3) / amu


    def calculate_ex_energy_grid_tensor(self,eigenvalues1,eigenvalues2,couplings1,couplings2,n_excits1,n_excits2,k_weights,n_points):



        dimension = np.shape(couplings1[0])[1]
        n_k_points = len(couplings1)


        ex_energy_grid = np.linspace(0,self.e_cutoff,n_points) #eV
        ex_energy_grid_tensor = np.zeros((dimension,dimension,n_points))
        ex_energy_grid_dx = ex_energy_grid[1]-ex_energy_grid[0]
        

        for i_spin in range(0,self.n_spin):

            for i_k_point in range (n_k_points):

                spin_k_factor = k_weights[i_k_point] * 2/self.n_spin

                eigenvalues1_k = eigenvalues1[i_k_point]
                couplings1_k = couplings1[i_k_point]
                n_excits1_k = n_excits1[i_k_point]

                eigenvalues2_k = eigenvalues2[i_k_point]
                couplings2_k = couplings2[i_k_point]
                n_excits2_k = n_excits2[i_k_point]

                # We take minimum number of excitations to prevent OOB error but this is funky
                n_excits_min = np.min(np.array((n_excits1_k,n_excits2_k)))

                fermi_occupation_factors = np.zeros((n_excits1_k))

                for i_ex in range(n_excits_min): # Approximation: use fermi occupation for structure 1
                    # We could analyse the fermi occupation for other structure and take average? or compare size

                    fermi_occupation_factors[i_ex] = fermi_pop(eigenvalues1_k[i_ex,0],self.chem_pot,self.temp) - \
                                                     fermi_pop(eigenvalues1_k[i_ex,1],self.chem_pot,self.temp)

                for i_coord in range(dimension):

                    for j_coord in range(dimension):

                        # .... discretize excitations
                        # don't need a separate loop over excitations in each structure since we
                        # assume the KS state energies do not change much, not sure how valid this 
                        for i_ex in range(n_excits1_k):
                            
                            ex1_e = eigenvalues1_k[i_ex,1] - eigenvalues1_k[i_ex,0]
                            ex2_e = eigenvalues2_k[i_ex,1] - eigenvalues2_k[i_ex,0]

                            delta = gaussian_function(ex1_e,ex_energy_grid[:],self.sigma) * \
                                    gaussian_function(ex2_e,ex_energy_grid[:],self.sigma)
                            norm = np.sum(delta)
                            norm = norm * ex_energy_grid_dx

                            if norm!=0.0:
                                ex_energy_grid_tensor[i_coord,j_coord,:] += fermi_occupation_factors[i_ex] * \
                                    np.real(np.conjugate(couplings1_k[i_ex,i_coord]) * \
                                    couplings2_k[i_ex,j_coord]) * \
                                    delta/norm * \
                                    spin_k_factor
  
        return ex_energy_grid, ex_energy_grid_tensor


    def evaluate_time_tensor(self,ex_energy_grid,ex_energy_grid_tensor,time_cutoff,n_time_points):

        dimension = np.shape(ex_energy_grid_tensor)[0]
        time_axis = np.linspace(0,time_cutoff,n_time_points)
        time_tensor = np.zeros((dimension,dimension,len(time_axis)))
        d_tau = time_axis[1]-time_axis[0]
        frequency_cutoff = self.e_cutoff/hbar

        # ..... Fourier transfoorm with cos


        for i_coord in range(dimension):

            for j_coord in range(dimension):

                for t,tau in enumerate(time_axis):
                    
                    integrand = np.zeros_like(ex_energy_grid)
                    
                    # Ignore e=0
                    integrand[1:] = ex_energy_grid_tensor[i_coord,j_coord,1:] * \
                        np.cos(ex_energy_grid[1:]/hbar * tau) / ex_energy_grid[1:]

                    time_tensor[i_coord,j_coord,t] =  simps(integrand,ex_energy_grid)



                    # ..... Convolute with sinc factor due to finite cutoff

                    sinc_factor = frequency_cutoff*np.sinc((tau*frequency_cutoff)/2)    

                    #exp_factor = np.cos(0.5*times_up*cutoff_freq) + (1j * np.sin(0.5*times_up*cutoff_freq))
                    exp_factor_real = np.cos(0.5*tau*frequency_cutoff)
    
                    convolute_factor = sinc_factor*exp_factor_real
                    
                    time_tensor[i_coord,j_coord,t] = np.convolve(time_tensor[i_coord,j_coord,t],convolute_factor,mode='same')*d_tau

        return time_axis, time_tensor


    def evaluate_time_tensor_no_convolute(self,ex_energy_grid,ex_energy_grid_tensor,time_cutoff,n_time_points):

        dimension = np.shape(ex_energy_grid_tensor)[0]
        time_axis = np.linspace(0,time_cutoff,n_time_points)
        time_tensor = np.zeros((dimension,dimension,len(time_axis)))
        d_tau = time_axis[1]-time_axis[0]
        frequency_cutoff = self.e_cutoff/hbar

        # ..... Fourier transfoorm with cos


        for i_coord in range(dimension):

            for j_coord in range(dimension):

                for t,tau in enumerate(time_axis):
                    
                    integrand = np.zeros_like(ex_energy_grid)
                    
                    # Ignore e=0
                    integrand[1:] = ex_energy_grid_tensor[i_coord,j_coord,1:] * \
                        np.cos(ex_energy_grid[1:]/hbar * tau) / ex_energy_grid[1:]

                    time_tensor[i_coord,j_coord,t] =  simps(integrand,ex_energy_grid)


class friction_gamma_parser_2022():
    def __init__(self,aims_file,gamma_files):
        self.chem_pot = self.parse_chem_pot(aims_file)
        self.chem_pot_response = self.parse_chem_pot_response(aims_file)
        #gamma_files.sort()
        max_k,gamma_files_sorted = self.sort_gamma_files(gamma_files)
        self.friction_masses = self.parse_fiction_masses(aims_file)
        coords,eis,ejs,couplings,kweights = self.parse_gamma_couplings(gamma_files_sorted,len(self.friction_masses))

        #self.ks = np.array(ks)
        self.max_k = max_k
        self.coords = np.array(coords)
        self.eis = eis
        self.ejs = ejs
        self.couplings = couplings
        self.kweights = np.array(kweights)



    def sort_gamma_files(self,gamma_files):

        gamma_files_sorted = np.array(gamma_files)

        k_no_list = []
        for i, gf in enumerate(gamma_files):

            if "spin" in gf:
                k_no = int(gf.replace('friction_gamma_k','').replace('.out','').replace("_spin_001",""))
            else:
                k_no = int(gf.replace('friction_gamma_k','').replace('.out',''))
            k_no_list.append(k_no)

        idx = np.argsort(k_no_list)

        gamma_files_sorted = gamma_files_sorted[idx]

        max_k = np.max(k_no_list)


        return max_k,gamma_files_sorted

    def parse_chem_pot(self,aims_file):
        chem_pot = 0
        with open(aims_file, "r") as af:
            for line in af:
                if '**FRICTION**' in line:
                    break
                if '| Chemical potential (Fermi level):' in line:
                    chem_pot = float(line.split()[-2])
        return chem_pot

    def parse_chem_pot_response(self,aims_file):
        chem_pot_response = []
        with open(aims_file, "r") as af:
            for line in af:
                if ' | Response of chemical potential (Fermi level):' in line:
                    chem_pot_response.append(float(line.split()[-2]))

        chem_pot_response = np.array(chem_pot_response)
        
        return chem_pot_response

    
    def parse_fiction_masses(self,aims_file):
        friction_atoms = []
        with open(aims_file, "r") as af:
            read = False
            for line in af:
                if 'The contents of geometry.in will be repeated verbatim below' in line:
                    read=True
                if "moment" in line:
                    continue
                if 'calculate_friction .true.' in line:
                    if (len(element)>2):
                        continue
                    else:
                        friction_atoms.append(element)
                if read:
                    try:
                        element = line.split()[-1]
                    except:
                        continue
        
        a = Atoms(symbols=friction_atoms)
        friction_masses = a.get_masses()

        return friction_masses

    def parse_gamma_couplings(self,gamma_files,n_f_atoms):
        #Parsing a limited amount to reduce memory usage
        print("--- %s Start parser ---" % (time.time() - start_time))
        ks = [] 
        coords = []
        # eis = [[0 for columns in range(0)] for i_rows in range(len(gamma_files))]
        eis = [[] for i_rows in range(len(gamma_files))]

        ejs = [[] for i_rows in range(len(gamma_files))]
        #ejs = [[0 for columns in range(0)] for i_rows in range(len(gamma_files))]
        couplings = [[0 for columns in range(0)] for i_rows in range(len(gamma_files))]
        kweights = []
        for i,file in enumerate(gamma_files):
            #print(file)
            couplings_k = []
            eis_k = []
            ejs_k = []

            with open(file, "r") as f:
                read_gamma = False

                for line in f:
                    if '| k-point:' in line:
                        #k = int(line.split()[2])
                        kweight = float(line.split()[-1])
                        kweights.append(kweight)
                    elif 'Friction component ' in line:
                        read_gamma = True
                        coord = int(line.split()[-1])-1
                        coords.append(coord)
                    elif read_gamma:
                        
                        ei = float(line.split()[0])

                        # if ei > 0.:
                        #     continue

                        ej = float(line.split()[1])
                        
                        # Trim eigenvalues
                        if ei>ej:
                            continue

                        #if ei == ej:
                        #     # if abs(ei-self.chem_pot) > 0.05:
                        #     if not (ei==self.chem_pot):
                        #    continue


                        # if (ej-ei)>1.0:
                        #     continue

                        # if (ej-ei)<1e-30:
                        #     continue

                        # if ei!=ej:
                        #     continue

                        # if ei > 2.8 or ej > 2.8:
                        #     continue

                        # if ei < -9.5 or ej < -9.5:
                        #     continue
                        # if ei == ej:
                        #     continue
                        # if abs(ei-ej)<1e-30:
                        #     continue

                        # if (abs(ei-ej)<1e-2 and ei!=ej):
                        #      continue
                        # if ej-ei > 3.0:
                        #     continue

                        # if ej < -6:
                        #     continue


                        coupling = complex(float(line.split()[2]),float(line.split()[3]))

                             
                        couplings_k.append(coupling)

                        if coord==1: #Eigenvalues same for each coomponent
                            eis_k.append(ei)
                            ejs_k.append(ej)

                couplings[i] = couplings_k
                eis[i] = eis_k
                ejs[i] = ejs_k

                        
        print("--- %s End parser ---" % (time.time() - start_time))
        return coords,eis,ejs,couplings,kweights

class friction_tensor_2022():
    #So ks are base zero (i.e first k_point = 0)
    #ks is just a list of what k point each excitation is associated with
    def __init__(self,parser,temp,sigma,nspin):
        self.couplings = parser.couplings
        self.coords = parser.coords
        self.eis = parser.eis
        self.ejs = parser.ejs
        #self.ks = parser.ks
        self.max_k = parser.max_k
        self.chem_pot = parser.chem_pot
        self.temp = temp
        self.kweights=parser.kweights
        self.sigma = sigma
        self.friction_masses = parser.friction_masses
        self.nspin = nspin

    def calc_tensor(self,mode='default',order=2,perturbing_energy=0.,fermi_correction='none',chemical_potential=None):
        print("--- %s Start calc tensor ---" % (time.time() - start_time))
        #ks = self.ks
        coords = self.coords
        ndim = np.max(coords)+1
        #max_k = np.max(ks)
        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings
        tensor = np.zeros((ndim,ndim))

        if chemical_potential is None:
            chem_pot = self.chem_pot
        else:
            chem_pot = chemical_potential


        temp = self.temp
        nspin = self.nspin

        if perturbing_energy > 0:
            if mode in ['default','double_delta','prb_print','double_delta_methfessel_paxton']:
                raise ValueError('mode can only deal with 0 perturbing energy')
            

        #Loop over k_points 
        for i_k,k in enumerate(range(1,self.max_k+1)):

            kw = self.kweights[k-1]
            es = np.array(ejs[i_k])-np.array(eis[i_k])

            if mode == 'no_norm':
                gaussian_delta_factor = ((gaussian_function(es,perturbing_energy,self.sigma))/(perturbing_energy))*kw*2/nspin 
            else:
                gaussian_delta_factor = ((gaussian_function(es,perturbing_energy,self.sigma))/(es))*kw*2/nspin 
            fermi_factor = fermi_pop(np.array(eis[i_k]),chem_pot,temp)-fermi_pop(np.array(ejs[i_k]),chem_pot,temp)

            couplings_k = np.array_split(np.array(couplings[i_k]), ndim)

            e_f_changes = [0,0,-4.22975078227829,0,0,5.86890516634832]
            e_f_changes = np.array(e_f_changes)

            omegas = es/hbar

            for i in range(ndim):
                #i_idx = np.where((coords == i) & (ks == k))[0]
                for j in range(ndim):
                    if j < i:
                        continue
                    #j_idx = np.where((coords == j) & (ks == k))[0]

                    nacs = np.conjugate(couplings_k[i])*couplings_k[j]

                    # for c_t,coupling in enumerate(nacs):
                    #     if eis[i_k][c_t] == eis[i_k][c_t]:
                    #         # nacs[c_t] = np.real(nacs[c_t]) + 0j
                    #         nacs[c_t] = np.conjugate((couplings_k[i][c_t] - e_f_changes[i])) * (couplings_k[j][c_t] - e_f_changes[j])


                    if mode=='default': # With additional normalisation described in RJM PRB 2016
                        # We no longer like this
                
                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor*gaussian_delta_factor/gaussian_norm(es,self.sigma))

                    elif mode=='double_delta':
              
                        tensor[i,j] += np.sum(nacs*\
                            (gaussian_function(np.array(eis[i_k])-chem_pot,0,self.sigma)*gaussian_function(np.array(ejs[i_k])-chem_pot,0,self.sigma))*\
                            kw*2/nspin)




                    elif mode=='occ_delta_only':
              
                        tensor[i,j] += np.sum(
                            (gaussian_function(np.array(eis[i_k])-chem_pot,0,self.sigma))*\
                            kw*2/nspin)


                    elif mode=='unocc_delta_only':
            
                        tensor[i,j] += np.sum(gaussian_function(np.array(ejs[i_k])-chem_pot,0,self.sigma)*\
                            kw*2/nspin)


                    elif mode=='double_delta_mauri':
              
                        tensor[i,j] += np.sum(nacs*\
                            (gaussian_function_mauri(np.array(eis[i_k]),0,self.sigma)*gaussian_function_mauri(np.array(ejs[i_k]),0,self.sigma))*\
                            kw*2/nspin)


                    elif mode=='novko':
              
                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor*(gaussian_derivative(np.array(eis[i_k]),np.array(ejs[i_k]),self.sigma))*\
                            kw*2/nspin)

                    elif mode=='double_delta_qe':

                            tensor[i,j] += np.sum(nacs*\
                            (gaussian_function_qe(eis[i_k],ejs[i_k],chem_pot,self.sigma))* \
                            kw*2/nspin)

                    elif mode=='prb_print':
            
                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor/(es)\
                            *(gaussian_function(es,0,self.sigma)/gaussian_norm2(es,self.sigma))*kw*2/nspin)


                    elif mode=='no_norm':

                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor*gaussian_delta_factor)

                    elif mode=='no_norm_approx':

                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor*gaussian_delta_factor)


                    elif mode=='fermi_derivative':
           
                        tensor[i,j] += np.sum(nacs*\
                            fermi_derivative(es,chem_pot,temp)\
                            *(gaussian_function(es,0,self.sigma))*kw*2/nspin)

                    elif mode == 'laplace':

                        omega_factors = (((perturbing_energy/hbar)*omegas)/(((perturbing_energy/hbar)**2)+omegas**2))
                        tensor[i,j] += np.sum((nacs/(es**2))*fermi_factor*omega_factors*kw*2/nspin)

                    elif mode=='gaussian2_no_norm':

                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor/(es)\
                            *(gaussian_function2(es,perturbing_energy,self.sigma))*kw*2/nspin)
                            
                    elif mode=='lorentzian': #no additional normalisation
            
                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor/(es)\
                            *(lorentzian_function(es,perturbing_energy,self.sigma))*kw*2/nspin)

                    elif mode=='methfessel_paxton': #no additional normalisation
                 
                        tensor[i,j] += np.sum(nacs*\
                            fermi_factor/(es)\
                            *(methfessel_paxton_function(es,perturbing_energy,self.sigma,order))*kw*2/nspin)

                    elif mode=='double_delta_methfessel_paxton':
                
                        tensor[i,j] += np.sum(nacs*\
                            (methfessel_paxton_function(np.array(eis[i_k])-chem_pot,0,self.sigma,order)*methfessel_paxton_function(np.array(ejs[i_k])-chem_pot,0,self.sigma,order))*\
                            kw*2/nspin)

                    else:
                        print('No viable tensor mode selected')
                    
        #tensor *= hbar*np.pi*2 
        #CLB 2021: don't think should have factor 2, removing it gives right tensor value
        tensor *= hbar*np.pi

        for i in range(ndim):
            mass_i = self.friction_masses[i // 3]
            for j in range(ndim):
                if j < i:
                    continue
                mass_j = self.friction_masses[j // 3]
                tensor[i,j] = tensor[i,j]/np.sqrt(mass_i*mass_j)
                tensor[j,i] = tensor[i,j]

        print("--- %s End calc tensor ---" % (time.time() - start_time))
        return tensor*ps
# couplings  [2,n_spin,n_k_points,dimension,dimension,3,max_n_excits]

# gamma_files = ['/Users/u1865573/Downloads/friction_gamma_k001.out','/Users/u1865573/Downloads/friction_gamma_k002.out']

# a = friction_output_parser_2021()
# a.parse_gamma_couplings(gamma_files)
