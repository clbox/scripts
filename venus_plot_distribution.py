import numpy as np
from ase.units import _hbar, J, s, fs
from ase import Atoms
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys

filenames = sys.argv[1:]

nstates = len(filenames)

fig, ax = plt.subplots(nstates, 1, sharex='all')#,sharey='all')
ymax=2.5

empty_ones = []
ntraj_list = []
state_list=[]
for i,filename in enumerate(filenames):

    print(filename)

    ntrajs = 0
    misc = []
    abs_e = []
    per_e = []

    with open(filename) as f:
        for line in f:

            if 'Trajectory' in line:
                ntrajs += 1
            
            elif 'Lifetime' in line:
                numbers = line.replace(',','')
                lifetime = float(numbers.split()[2])
                scat = float(numbers.split()[7])
                jf = int(numbers.split()[-1])

                misc.append([lifetime,scat,jf])

            elif 'Total' in line:
                numbers = line.replace(',','')
                d = float(numbers.split()[7])
                theta = float(numbers.split()[10])
                phi = float(numbers.split()[13])
                X = float(numbers.split()[16])
                Y = float(numbers.split()[19])
                Z = float(numbers.split()[22])
                tot = float(numbers.split()[-1])

                abs_e.append([d,theta,phi,X,Y,Z,tot])

            elif '%' in line:
                numbers = line.replace(',','')
                d = float(numbers.split()[6])
                theta = float(numbers.split()[9])
                phi = float(numbers.split()[12])
                X = float(numbers.split()[15])
                Y = float(numbers.split()[18])
                Z = float(numbers.split()[21])

                per_e.append([d,theta,phi,X,Y,Z])

    misc = np.array(misc)
    abs_e = np.array(abs_e)
    per_e = np.array(per_e)

    if 'trapped' in filename:
        ax[0].boxplot(abs_e,showfliers=False)
        ax[0].text(s=filename,x=2,y=np.max(abs_e)-0.5*np.max(abs_e))
        ax[0].text(s='ntrajs = '+str(ntrajs),x=5,y=np.max(abs_e)-0.5*np.max(abs_e))
        ax[0].set_ylim(0,np.max(abs_e))
    else:
        ax[i+1].set_ylim(0,ymax)
        ax[i+1].boxplot(abs_e,showfliers=False)
        ax[i+1].text(s=filename,x=2,y=ymax-0.5*ymax)
        ax[i+1].text(s='ntrajs = '+str(ntrajs),x=5,y=ymax-0.5*ymax)
        print('Average lifetime / fs ' +str(np.average(misc[:,0])))
        print('Average scattering angle ' + str(np.average(misc[:,1])))
        print('Average final rotational state ' + str(np.average(misc[:,2])))
        ntraj_list.append(ntrajs)
        state_list.append(int(filename.split('.')[0]))

plt.xticks([1,2,3,4,5,6,7],['d',r'$\theta$',r'$\phi$','X','Y','Z','Total'])


fig.set_figheight(10*nstates/10)
fig.set_figwidth(7)
fig.text(0.5, 0.05, "Mode", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Energy loss / eV', va='center', rotation='vertical',fontsize=15)
fig.savefig('summary.pdf',transparent=True,bbox_inches='tight')



ntraj_list = np.array(ntraj_list)
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
ax.plot(state_list,ntraj_list/np.sum(ntraj_list),'.-')
fig.set_figheight(4)
fig.set_figwidth(7)
fig.text(0.5, 0.05, "Nf", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Probability', va='center', rotation='vertical',fontsize=15)
fig.savefig('probability.pdf',transparent=True,bbox_inches='tight')