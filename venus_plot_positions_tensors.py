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
from ase.io import read
from ase.io.trajectory import Trajectory


traj_no = sys.argv[1]

time_axis = np.loadtxt('timeaxis_{}.txt'.format(traj_no))
projected_tensors = vjp.read_array('projected_tensors{}.txt'.format(traj_no))
traj = Trajectory('trajectory_{}.traj'.format(traj_no))

# time_axis = np.loadtxt(sys.argv[1])
# projected_tensors = vjp.read_array(sys.argv[2])
# traj = Trajectory(sys.argv[3])

labels2 = [r'$\Lambda_{rr}$',r'$\Lambda_{\theta\theta}$',r'$\Lambda_{\phi\phi}$',r'$\Lambda_{XX}$',r'$\Lambda_{YY}$',r'$\Lambda_{ZZ}$']
labels1 = [r'$\dot{R}_{rr}$',r'$\dot{R}_{\theta\theta}$',r'$\dot{R}_{\phi\phi}$',r'$\dot{R}_{XX}$',r'$\dot{R}{YY}$',r'$\dot{R}_{ZZ}$']
line_args = {'linewidth' : 1.0}
colours1 = ['red','dodgerblue','black','magenta','orange','orange']
colours2 = ['purple','dodgerblue','black','magenta','orange','orange']

o_pos=[]
n_pos=[]

for i,ts in enumerate(traj):
    atoms = traj[i]
    positions = atoms.get_positions()
    #Z positions
    o_pos.append(positions[0,2])
    n_pos.append(positions[1,2])



ndim = np.shape(projected_tensors)[1]

fig, ax1 = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)
ax2 = ax1.twinx()
for i in range(ndim):
    if i in [1,2,3,4,5]:
        continue
    if i <=2:
        linestyle = '-'
    else:
        linestyle='-'
   
    ax2.plot(time_axis,projected_tensors[:,i,i],label=labels2[i],linestyle='-',**line_args,color=colours2[i])

ax1.plot(time_axis,o_pos,label='Oxygen',linestyle='-.',**line_args,color='red')
ax1.plot(time_axis,n_pos,label='Nitrogen',linestyle='--',**line_args,color='navy')

for ax in [ax1,ax2]:
    font='Arial'
    for tick in ax.get_xticklabels():
        tick.set_fontname(font)
    for tick in ax.get_yticklabels():
        tick.set_fontname(font)

    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))


ax1.yaxis.set_minor_locator(MultipleLocator(0.2))
ax1.yaxis.set_major_locator(MultipleLocator(1))

ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.xaxis.set_major_locator(MultipleLocator(0.1))

ax1.set_xlim(left=0.0,right=0.4)
ax1.set_ylim(bottom=1.0,top=7)


ax2.set_ylim(bottom=0.0,top=1.2)
ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
ax2.yaxis.set_major_locator(MultipleLocator(0.2))

fig.set_figheight(2.0)
fig.set_figwidth(3.25)

plt.gcf().subplots_adjust(left=0.3,bottom=0.3,right=0.8)


ax1.set_xlabel('Time / ps',fontsize=12,fontname=font)#,color='white')
ax1.set_ylabel(r'$R_z$ / $\mathrm{ \AA{}}$',fontsize=12,fontname=font,color='black')#,color='white')
ax2.set_ylabel(r'$\Lambda_{ij}$ / ps$^{-1}$',fontsize=12,fontname=font,color='purple')#,color='white')
fig.savefig('pos_tensors_{}.pdf'.format(traj_no),transparent=True)#,bbox_inches='tight')


fig.legend(ncol=6,fontsize=12,fancybox=True,framealpha=0)
fig.set_figwidth(7)
ax1.remove()
ax2.remove()
fig.savefig('legend.pdf',transparent=True)#,bbox_inches='tight')