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

def gaussian_norm(x0, s):
    return 0.5 * (1-special.erf((-x0/s)*(1/np.sqrt(2))))

def gaussian_norm2(x0, s):
    return 0.5 * (1-special.erf((-x0)))

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

class friction_gamma_parser():
    def __init__(self,aims_file,gamma_files):
        self.chem_pot = self.parse_chem_pot(aims_file)
        #gamma_files.sort()
        ks,coords,eis,ejs,couplings,kweights = self.parse_gamma_couplings(gamma_files)
        self.ks = np.array(ks)
        self.coords = np.array(coords)
        self.eis = np.array(eis)
        self.ejs = np.array(ejs)
        self.couplings = np.array(couplings)
        self.kweights = np.array(kweights)
        self.friction_masses = self.parse_fiction_masses(aims_file)

    def parse_chem_pot(self,aims_file):
        chem_pot = 0
        with open(aims_file, "r") as af:
            for line in af:
                if '**FRICTION**' in line:
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
        ks = [] 
        coords = []
        eis = []
        ejs = []
        couplings = []
        kweights = []
        for file in gamma_files:
            print(file)
            with open(file, "r") as f:
                read_gamma = False
                for line in f:
                    if '| k-point:' in line:
                        k = int(line.split()[2])
                        kweight = float(line.split()[-1])
                        kweights.append(kweight)
                    elif 'Friction component ' in line:
                        read_gamma = True
                        coord = int(line.split()[-1])-1
                    elif read_gamma:
                        ei = float(line.split()[0])
                        ej = float(line.split()[1])
                        coupling = complex(float(line.split()[2]),float(line.split()[3]))
                        
                        if ei>ej:
                            continue
                        if ei == ej:
                            continue
                        if abs(ei-ej)<1e-30:
                            continue
                        if ej-ei > 3.0:
                            continue
                        if ei > 0.:
                            continue
                        if ej < -6:
                            continue
                        #ks.append(k-1)
                        ks.append(k)
                        coords.append(coord)
                        eis.append(ei)
                        ejs.append(ej)
                        couplings.append(coupling)
        print("--- %s End parser ---" % (time.time() - start_time))
        return ks,coords,eis,ejs,couplings,kweights

class friction_tensor():
    #So ks are base zero (i.e first k_point = 0)
    #ks is just a list of what k point each excitation is associated with
    def __init__(self,parser,temp,sigma,nspin):
        self.couplings = parser.couplings
        self.coords = parser.coords
        self.eis = parser.eis
        self.ejs = parser.ejs
        self.ks = parser.ks
        self.chem_pot = parser.chem_pot
        self.temp = temp
        self.kweights=parser.kweights
        self.sigma = sigma
        self.friction_masses = parser.friction_masses
        self.nspin = nspin

    def calc_tensor(self,mode='default',order=2,perturbing_energy=0.):
        print("--- %s Start calc tensor ---" % (time.time() - start_time))
        ks = self.ks
        coords = self.coords
        ndim = np.max(coords)
        max_k = np.max(ks)
        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings
        tensor = np.zeros((ndim+1,ndim+1))
        chem_pot = self.chem_pot
        temp = self.temp
        nspin = self.nspin
        for k in range(1,max_k+1):
            #Loop over k_points 
            #kw = self.kweights[k]
            #i_kw = np.where(ks == k)[0]
            kw = self.kweights[k-1]
            
            for i in range(ndim+1):
                i_idx = np.where((coords == i) & (ks == k))[0]
                for j in range(ndim+1):
                    if j < i:
                        continue
                    j_idx = np.where((coords == j) & (ks == k))[0]

                    es = ejs[j_idx]-eis[i_idx]

                    if mode=='default':
                        if perturbing_energy!=0.:
                            print('default mode can only deal with 0 perturbing energy')
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                        (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                            *(gaussian_function(es,0,self.sigma)/gaussian_norm(es,self.sigma))*kw*2/nspin)

                    elif mode=='double_delta':
                        if perturbing_energy!=0.:
                            print('double delta mode can only deal with 0 perturbing energy')
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                            (gaussian_function(eis[i_idx]-chem_pot,0,self.sigma)*gaussian_function(ejs[j_idx]-chem_pot,0,self.sigma))*\
                            kw*2/nspin)

                    # elif mode=='double_delta_half_sigma':
                    #     if perturbing_energy!=0.:
                    #         print('double delta mode can only deal with 0 perturbing energy')
                    #     tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                    #         (gaussian_function(eis[i_idx]-chem_pot,0,self.sigma/2)*gaussian_function(ejs[j_idx]-chem_pot,0,self.sigma/2))*\
                    #         kw*2/nspin)

                    elif mode=='prb_print':
                        if perturbing_energy!=0.:
                            print('prb print mode can only deal with 0 perturbing energy')
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                        (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                            *(gaussian_function(es,0,self.sigma)/gaussian_norm2(es,self.sigma))*kw*2/nspin)

                    elif mode=='no_norm':
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                        (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                            *(gaussian_function(es,perturbing_energy,self.sigma))*kw*2/nspin)

                    elif mode=='gaussian2_no_norm':
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                        (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                            *(gaussian_function2(es,perturbing_energy,self.sigma))*kw*2/nspin)
                            
                    elif mode=='lorentzian': #no additional normalisation
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                        (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                            *(lorentzian_function(es,perturbing_energy,self.sigma))*kw*2/nspin)

                    elif mode=='methfessel_paxton': #no additional normalisation
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                        (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                            *(methfessel_paxton_function(es,perturbing_energy,self.sigma,order))*kw*2/nspin)

                    elif mode=='double_delta_methfessel_paxton':
                        if perturbing_energy!=0.:
                            print('double delta mode can only deal with 0 perturbing energy')
                        tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                            (methfessel_paxton_function(eis[i_idx]-chem_pot,0,self.sigma,order)*methfessel_paxton_function(ejs[j_idx]-chem_pot,0,self.sigma,order))*\
                            kw*2/nspin)
                    else:
                        print('No viable tensor mode selected')
                    
        #tensor *= hbar*np.pi*2 
        #CLB 2021: don't think should have factor 2, removing it gives right tensor value
        tensor *= hbar*np.pi

        for i in range(ndim+1):
            mass_i = self.friction_masses[i // 3]
            for j in range(ndim+1):
                if j < i:
                    continue
                mass_j = self.friction_masses[j // 3]
                tensor[i,j] = tensor[i,j]/np.sqrt(mass_i*mass_j)
                tensor[j,i] = tensor[i,j]

        print("--- %s End calc tensor ---" % (time.time() - start_time))
        return tensor*ps

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
                        # assumption is but I think its already implicit within MDEF anyway.
                        for i_ex in range(n_excits_min):
                            
                            ex1_e = eigenvalues1_k[i_ex,1] - eigenvalues1_k[i_ex,0]
                            ex2_e = eigenvalues2_k[i_ex,1] - eigenvalues2_k[i_ex,0]

                            ex_energy_grid_tensor[i_coord,j_coord,:] += fermi_occupation_factors[i_ex] * \
                                np.real(np.conjugate(couplings1_k[i_ex,i_coord]) * \
                                couplings2_k[i_ex,j_coord]) * \
                                gaussian_function(ex1_e,ex_energy_grid[:],self.sigma) * \
                                gaussian_function(ex2_e,ex_energy_grid[:],self.sigma) * \
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



# couplings  [2,n_spin,n_k_points,dimension,dimension,3,max_n_excits]

# gamma_files = ['/Users/u1865573/Downloads/friction_gamma_k001.out','/Users/u1865573/Downloads/friction_gamma_k002.out']

# a = friction_output_parser_2021()
# a.parse_gamma_couplings(gamma_files)
