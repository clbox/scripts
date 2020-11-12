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
    fermi_pop = (1./(np.exp((x-x0)/(boltzmann_kB*T))+1))
    return fermi_pop

def gaussian_function(x,x0,s):
    return 1./np.sqrt(2*np.pi*s*s)*np.exp(-0.5*(((x-x0) / (s))**2))

def gaussian_norm(x0, s):
    return 0.5 * (1-special.erf((-x0/s)*(1/np.sqrt(2))))

class friction_gamma_parser():
    def __init__(self,aims_file,gamma_files):
        self.chem_pot = self.parse_chem_pot(aims_file)
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
            with open(file, "r") as f:
                read_gamma = False
                for line in f:
                    if '| k-point:' in line:
                        k = int(line.split()[2])
                        kweight = float(line.split()[-1])
                    elif 'Friction component ' in line:
                        read_gamma = True
                        coord = int(line.split()[-1])-1
                    elif read_gamma:
                        ei = float(line.split()[0])
                        ej = float(line.split()[1])
                        coupling = complex(float(line.split()[2]),float(line.split()[3]))
                        
                        ks.append(k-1)
                        coords.append(coord)
                        eis.append(ei)
                        ejs.append(ej)
                        couplings.append(coupling)
                        kweights.append(kweight)
        print("--- %s End parser ---" % (time.time() - start_time))
        return ks,coords,eis,ejs,couplings,kweights

class friction_tensor():
    def __init__(self,ks,coords,eis,ejs,couplings,kweights,chem_pot,friction_masses,temp,sigma,nspin):
        self.couplings = couplings
        self.coords = coords
        self.eis = eis
        self.ejs = ejs
        self.ks = ks
        self.chem_pot = chem_pot
        self.temp = temp
        self.kweights=kweights
        self.sigma = sigma
        self.friction_masses = friction_masses
        self.nspin = nspin

    def calc_tensor(self):
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

        for k in range(max_k):
            kw = self.kweights[k]
            for i in range(ndim+1):
                i_idx = np.where((coords == i) & (ks == k))[0]
                for j in range(ndim+1):
                    if j < i:
                        continue
                    j_idx = np.where((coords == j) & (ks == k))[0]
                    es = ejs[j_idx]-eis[i_idx]
                    tensor[i,j] += np.sum(np.conjugate(couplings[i_idx])*couplings[j_idx]*\
                    (fermi_pop(eis[i_idx],chem_pot,temp)-fermi_pop(ejs[j_idx],chem_pot,temp))/(es)\
                        *(gaussian_function(es,0,self.sigma)/gaussian_norm(es,self.sigma))*kw*2/nspin)
        tensor *= hbar*np.pi*2

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


class calc_time_tensor():
    #discretize ei and ejs with gaussian functions of smearing width sigma
    def __init__(self,ks,coords,eis,ejs,couplings,kweights,chem_pot,nspin,min_e,max_e,npoints,friction_masses,temp=300,sigma=0.01):
        self.couplings = couplings
        self.coords = coords
        self.eis = eis
        self.ejs = ejs
        self.ks = ks
        self.chem_pot = chem_pot
        self.kweights=kweights
        self.temp = temp
        self.sigma = sigma
        self.nspin = nspin
        self.min_e = min_e
        self.max_e = max_e
        self.npoints = npoints
        self.ndim = np.max(coords)
        self.friction_masses = friction_masses

    def get_e_grid(self):
        min_e = self.min_e
        max_e = self.max_e
        npoints = self.npoints
        e_grid =  np.linspace(min_e,max_e,npoints)
        return e_grid

    def get_h_grid(self):
        min_e = self.min_e
        max_e = self.max_e
        npoints = self.npoints
        h_grid = np.linspace(min_e,max_e,npoints)
        return h_grid
    
    def get_ex_grid(self):
        min_e = self.min_e
        max_e = self.max_e
        npoints = self.npoints
        h_grid = np.linspace(min_e,max_e,npoints)
        return h_grid

    def calc_A(self):
        print("--- %s Start calc A ---" % (time.time() - start_time))

        coords = self.coords
        ndim = self.ndim
        temp = self.temp
        chem_pot = self.chem_pot
        e_grid = self.get_e_grid()
        h_grid = self.get_h_grid()
        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings
        sigma = self.sigma

        A = np.zeros((ndim+1,ndim+1,len(e_grid),len(h_grid)),dtype=np.complex128)

        for n1 in range(ndim+1):
            for n2 in range(ndim+1):
                if n2< n1:
                    continue
                idx = np.where((coords == n1))[0]
                idx2 = np.where((coords == n2))[0]
                cs = couplings[idx]
                cs2 = couplings[idx2]
                eis_n = eis[idx]
                ejs_n = ejs[idx2]
                for i,c2 in enumerate(cs2):
                    c = np.conjugate(cs[i])
                    gauss_mat = np.outer(gaussian_function(eis_n[i],e_grid,sigma),gaussian_function(ejs_n[i],h_grid,sigma))\
                        *(fermi_pop(eis_n[i],chem_pot,temp)-fermi_pop(ejs_n[i],chem_pot,temp))
                    #print(np.max(gaussian_function(ejs_n[i],h_grid,0.01)))
                    A[n1,n2,:,:] +=gauss_mat*c*c2
        print("--- %s End calc A ---" % (time.time() - start_time))
        return A
    
    def calc_A2(self):
        print("--- %s Start calc A2 ---" % (time.time() - start_time))
        coords = self.coords
        ndim = self.ndim
        temp = self.temp
        chem_pot = self.chem_pot
        ex_grid = self.get_ex_grid()
        eis = self.eis
        ejs = self.ejs
        couplings = self.couplings
        sigma = self.sigma

        A = np.zeros((ndim+1,ndim+1,len(ex_grid)),dtype=np.complex128)

        for n1 in range(ndim+1):
            for n2 in range(ndim+1):
                if n2< n1:
                    continue
                idx = np.where((coords == n1))[0]
                idx2 = np.where((coords == n2))[0]
                cs = couplings[idx]
                cs2 = couplings[idx2]
                eis_n = eis[idx]
                ejs_n = ejs[idx2]
                for i,c2 in enumerate(cs2):
                    c = np.conjugate(cs[i])
                    gauss_mat = gaussian_function(ejs_n[i]-eis_n[i],ex_grid,sigma)*gaussian_function(ejs_n[i]-eis_n[i],ex_grid,sigma)\
                        /
                        #*(fermi_pop(eis_n[i],chem_pot,temp)-fermi_pop(ejs_n[i],chem_pot,temp))
                    A[n1,n2,:] +=gauss_mat*c*c2
        print("--- %s End calc A2 ---" % (time.time() - start_time))
        return A

    def mass_weight(self,tensor):
        ndim = self.ndim
        for i in range(ndim+1):
            mass_i = self.friction_masses[i // 3]
            for j in range(ndim+1):
                if j < i:
                    continue
                mass_j = self.friction_masses[j // 3]
                tensor[i,j] = tensor[i,j]/np.sqrt(mass_i*mass_j)
                tensor[j,i] = tensor[i,j]
        return tensor

    def evaluate_tensor(self):
        print("--- %s Start eval tensor ---" % (time.time() - start_time))
        ndim = self.ndim
        h_grid = self.get_h_grid()
        e_grid =  self.get_e_grid()
        A = self.calc_A()
        for i in range(len(e_grid)):
            for j in range(len(h_grid)):
                if (h_grid[j]-e_grid[i])==0:
                    A[:,:,i,j]=0
                else:
                    A[:,:,i,j] = A[:,:,i,j]/(h_grid[j]-e_grid[i])


        tensor = np.zeros((ndim+1,ndim+1))
        for n1 in range(ndim+1):
            for n2 in range(ndim+1):
                if n2< n1:
                    continue
                b = simps(A[n1,n2,:,:],e_grid,axis=0)
                tensor[n1,n2] = simps(b,h_grid)
                tensor[n2,n1] = tensor[n1,n2]

        tensor = self.mass_weight(tensor)
        tensor *= 2*ps*ps
        print("--- %s End eval tensor ---" % (time.time() - start_time))
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

                


