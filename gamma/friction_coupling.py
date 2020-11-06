import numpy as np
import sys
import os
from ase import Atoms
from scipy import special
from ase.units import _hbar, J, s, fs
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

        return ks,coords,eis,ejs,couplings,kweights

class friction_tensor():
    def __init__(self,ks,coords,eis,ejs,couplings,kweights,chem_pot,friction_masses,temp,sigma,nspin):
        self.weighted_couplings = couplings
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
        ks = self.ks
        coords = self.coords
        ndim = np.max(coords)
        max_k = np.max(ks)
        eis = self.eis
        ejs = self.ejs
        couplings = self.weighted_couplings
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


        return tensor*ps

                


