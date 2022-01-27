import numpy as np
import sys
import os
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)

import venus_traj_project as vjp
plt.style.use('clb_publication')


time_axis = np.loadtxt(sys.argv[1])

projected_tensors = vjp.read_array(sys.argv[2])
projected_velocities = vjp.read_array(sys.argv[3])

labels2 = [r'$\Lambda_{rr}$',r'$\Lambda_{\theta\theta}$',r'$\Lambda_{\phi\phi}$',r'$\Lambda_{XX}$',r'$\Lambda_{YY}$',r'$\Lambda_{ZZ}$']
labels1 = [r'$\dot{R}_{rr}$',r'$\dot{R}_{\theta\theta}$',r'$\dot{R}_{\phi\phi}$',r'$\dot{R}_{XX}$',r'$\dot{R}{YY}$',r'$\dot{R}_{ZZ}$']



ndim = np.shape(projected_velocities)[1]


projected_velocities = np.square(projected_velocities)

fig, ax1 = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)
ax2 = ax1.twinx()
for i in range(ndim):
    if i in [1,2,3,4,5]:
        continue
    if i <=2:
        linestyle = '-'
    else:
        linestyle='-'
    ax1.plot(time_axis,projected_velocities[:,i],label=labels1[i],marker='None',linewidth=1,linestyle='-',color='#1F90FF')
    ax2.plot(time_axis,projected_tensors[:,i,i],label=labels2[i],marker='None',linewidth=2,linestyle='--')

for ax in [ax1,ax2]:
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))


ax1.yaxis.set_minor_locator(MultipleLocator(0.0005))
ax1.yaxis.set_major_locator(MultipleLocator(0.001))

ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.xaxis.set_major_locator(MultipleLocator(0.1))

ax1.set_xlim(left=0.0,right=0.4)
ax1.set_ylim(bottom=0.0,top=0.001)


ax2.set_ylim(bottom=0.0,top=0.8)
ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
ax2.yaxis.set_major_locator(MultipleLocator(0.2))

fig.set_figheight(2.0)
fig.set_figwidth(3.25)

plt.gcf().subplots_adjust(left=0.3,bottom=0.3,right=0.8)


ax1.set_xlabel('Time / ps')#,color='white')
ax1.set_ylabel(r'$\mid\dot{R_{rr}}\mid ^{2}$ / $\mathrm{ \AA{}}^{2}$ fs$^{-2}$',color='#1F90FF')#,color='white')
ax2.set_ylabel(r'$\Lambda_{rr}$ / ps$^{-1}$',color='#FF4500')#,color='white')
fig.savefig('both_new.pdf',transparent=True)#,bbox_inches='tight')


fig.legend(ncol=6)
fig.set_figwidth(7)
ax1.remove()
ax2.remove()
fig.savefig('legend_new.pdf',transparent=True)#,bbox_inches='tight')