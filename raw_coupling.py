import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
import glob
from ase.io import read



output_dir = sys.argv[1]

nacs_files = (glob.glob(output_dir+'/*nacs-spectrum*.out'))

geo_file =  (glob.glob(output_dir+'/*geometry.in'))[0]
print('-------------'+geo_file+'----------------')
atoms = read(geo_file)

output = (glob.glob(output_dir+'/*aims.out'))[0]
print('-------------'+output+'----------------')


with open(output, "r") as f:
    for line in f:
        if '************************FRICTION**********************************' in line:
            break
        if '| Chemical potential (Fermi level):' in line:
            fermi = float(line.split()[-2])


results = {'e_occ' : [], 'e_unocc' : [], 'coupling' : [], 'k' : [], 'i_coord' : [], 'j_coord' : []}
nacs_files.sort()
print(nacs_files)
header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction","Orb"] #skip lines
for k,nacs in enumerate(nacs_files):
    
    with open(nacs, "r") as f:
        for line in f:

            if "Friction" in line:
                i_coord = int(line.split()[3])
                j_coord = int(line.split()[4])
                continue

            if any(x in line for x in header):
                continue

            results['e_occ'].append(float(line.split()[0]))
            results['e_unocc'].append(float(line.split()[1]))
            results['coupling'].append(float(line.split()[2]))
            results['k'].append(k)
            results['i_coord'].append(i_coord)
            results['j_coord'].append(j_coord)


e_occs = np.array(results['e_occ'])
e_unoccs = np.array(results['e_unocc'])
couplings = np.array(results['coupling'])
ks = np.array(results['k'])
i_coords = np.array(results['i_coord'])
j_coords = np.array(results['j_coord'])

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')
dimension = np.max(i_coords)

idx = np.argwhere((i_coords==1) & (j_coords == 1)).flatten()

# ax.scatter(e_occs[idx]-fermi,couplings[idx],marker='.')
# ax.scatter(e_unoccs[idx]-fermi,couplings[idx],marker='x')
excit_es = e_unoccs[idx]-e_occs[idx]
#ax.scatter(excit_es,couplings[idx])
bins = np.linspace(0,1.0,400)

# ax.hist(e_unoccs[idx]-fermi,color='maroon',bins=bins,label='Unoccupied')
# ax.hist(e_occs[idx]-fermi,color='dodgerblue',bins=bins,label='Occupied')



#ax.legend()

ax.hist(excit_es,bins=bins,color='navy')

ax.set_xlim(0,1.0)
ax.set_ylim(bottom=0)
#ax.set_xlabel("Excitation energy / eV")
ax.set_xlabel(r"$\epsilon_i - \epsilon_F$ / eV")
ax.set_ylabel(r'Frequency of excitations')
fig.set_figheight(3.0)
fig.set_figwidth(3.25)
fig.savefig('fig.pdf',transparent=True,bbox_inches='tight')