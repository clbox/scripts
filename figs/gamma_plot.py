import numpy as np
from scipy.interpolate import interp1d
from ase.units import _hbar, J, s, fs
from ase import Atoms
import math
import sys
hbar = _hbar * J * s 
ps = fs*1000

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import memory_eta_analysis as mea


#masses = mea.get_friction_masses('geometry.in',[16,17])

paths = sys.argv[1:]

labels = [r'C$_x$',r'C$_y$',r'C$_z$',r'O$_x$',r'O$_y$',r'O$_z$']

colours = ['red','blue','purple','green','yellow','black','cyan','orange']

bins,re_memory_kernel,im_memory_kernel,dimension,max_e = mea.read_memory_kernel(paths[0],treat_complex=False)
fig, ax = plt.subplots(dimension,1,sharex='all',sharey='all')

xmax = 1
ymax = np.max(re_memory_kernel[:,bins<xmax])
ymin = np.min(re_memory_kernel[:,bins<xmax])

for p,path in enumerate(paths):
    bins,re_memory_kernel,im_memory_kernel,dimension,max_e = mea.read_memory_kernel(path,treat_complex=False)
    for i in range(dimension):
        ax[i].plot(bins,re_memory_kernel[i,:],'-',color=colours[p],label='Re '+path.split('/')[0])
        #ax[i].plot(bins,im_memory_kernel[i,:],'--',color=colours[p],label='Im '+path.split('/')[0])
        ax[i].text(s=labels[i],x=0.1,y=0.55*ymax,fontsize=15)



ax[0].set_xlim(0,xmax)
ax[0].set_ylim(ymin,ymax)
ax[0].legend()
fig.set_figheight(7)
fig.set_figwidth(4)
fig.text(0.5, 0.01, r"Excitation energy / eV", ha='center',fontsize=15)
fig.text(-0.1, 0.5, r"$\sum_{\mathbf{k},v,v'} g_{\mathbf{k},v,v'} \tilde{\delta}(\epsilon_{v'}-\epsilon_{v}-\epsilon)$ / $\AA^{-1}$", va='center', rotation='vertical',fontsize=15)
fig.savefig('gamma.pdf',bbox_inches = "tight")