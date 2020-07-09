import numpy as np
import glob
import sys
import os
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

dir1 = sys.argv[1]

energy_files = ['./none/'+str(dir1)+'/energies.out','./collinear/'+str(dir1)+'/energies.out']

linestyles = ['-',':']
markers = ['o','^']
labels = [r'$E(S = 0)$',r'$E(S = 1)$']

dirs = ['none','collinear']

height_dirs = np.arange(0,15,1)

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')

with open(energy_files[1]) as f:
    lines = f.readlines()
    y = [float(line.split()[1]) for line in lines]
    y = np.array(y)
    lowest_energy = y[-1]

for i,filename in enumerate(energy_files):
    with open(filename) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]
        x = np.array(x)
        if sys.argv[1]=='ndown' or sys.argv[1]=='odown':
            x = x +2
        x = (x*0.5) + 0.1
        print(x)
        y = np.array(y)
        #normalise_y = y[-1]
        #y = y-normalise_y
        y = y - lowest_energy
        ax.plot(x,y,linestyle=linestyles[i],label=labels[i],color='gray',linewidth=2,marker=markers[i])

#PBE
ax.set_ylim( -19287077.6, -19287076.7)
#BEEF
#ax.set_ylim(-19298636,-19298630)
#NORMALISE
ax.set_ylim(-0.1,1.3)


ax2 = ax.twinx()
for i,dir in enumerate(dirs):

    heights = []
    tensors = []
    
    spin = str(i)

    for ii,height_dir in enumerate(height_dirs):

        try:
            #ft = np.loadtxt(dir+'/be/{:02d}/friction_tensor.out'.format(height_dir))
            ft = np.loadtxt(dir+'/{}/{:02d}/friction_tensor.out'.format(dir1,height_dir))
        except:
            continue
        
        heights.append((height_dir*0.5)+0.1)
        tensors.append(ft)


    heights = np.array(heights)
    tensors = np.array(tensors)
    
    #ax2.plot(heights,tensors[:,0,0],linestyle=linestyles[i],label=r'$\Lambda_{O_x O_x}(S=$'+spin+r'$)$',color='blue')
    #ax2.plot(heights,tensors[:,1,1],linestyle=linestyles[i],label=r'$\Lambda_{O_y O_y}(S=$'+spin+r'$)$',color='red')
    ax2.plot(heights,tensors[:,2,2],linestyle=linestyles[i],label=r'$\Lambda_{O_z O_z}(S=$'+spin+r'$)$',color='red',marker='.')



    #ax2.plot(heights,tensors[:,3,3],linestyle=linestyles[i],label=r'$\Lambda_{N_x N_x}(S=$'+spin+r'$)$',color='navy')
    #ax2.plot(heights,tensors[:,4,4],linestyle=linestyles[i],label=r'$\Lambda_{N_y N_y}(S=$'+spin+r'$)$',color='maroon')
    ax2.plot(heights,tensors[:,5,5],linestyle=linestyles[i],label=r'$\Lambda_{N_z N_z}(S=$'+spin+r'$)$',color='navy',marker='.')
    
    #ax2.plot(heights,tensors[:,1,1],linestyle=linestyles[i],label=r'$\Lambda_{\theta \theta}(S=$'+spin+r'$)$',color='orangered',marker='.')
    #ax2.plot(heights,tensors[:,2,2],linestyle=linestyles[i],label=r'$\Lambda_{\phi \phi}(S=$'+spin+r'$)$',color='palegreen',marker='.')
    #ax2.plot(heights,tensors[:,0,0],linestyle=linestyles[i],label=r'$\Lambda_{rr}(S=$'+spin+r'$)$',color='deepskyblue',marker='.')


ax2.set_ylim(0,1.6)
ax.set_xlim(0,7.1)

ax.legend(loc=2)
ax2.legend(loc=0,ncol=2)
fig.text(0.5, 0.01, r'COM height above the surface / $\AA$', ha='center',fontsize=15)
fig.text(0.01, 0.5, r"$\mathrm{E_{ads}}$ / eV", va='center',rotation='vertical',fontsize=15)
fig.text(1.0, 0.5, r"$\Lambda_{ij}$ / ps$^{-1}$", va='center',rotation='vertical',fontsize=15)
fig.set_figwidth(6)
fig.set_figheight(5)
fig.savefig('energies.pdf',transparent=True,bbox_inches='tight')







