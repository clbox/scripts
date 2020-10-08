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

import venus.venus_traj_project as vjp



labels2 = [r'$\Lambda_{rr}$',r'$\Lambda_{\theta\theta}$',r'$\Lambda_{\phi\phi}$',r'$\Lambda_{XX}$',r'$\Lambda_{YY}$',r'$\Lambda_{ZZ}$']
labels1 = [r'$\dot{R}_{rr}$',r'$\dot{R}_{\theta\theta}$',r'$\dot{R}_{\phi\phi}$',r'$\dot{R}_{XX}$',r'$\dot{R}{YY}$',r'$\dot{R}_{ZZ}$']
line_args = {'linewidth' : 0.8}
colours1 = ['red','dodgerblue','black','magenta','orange','orange']
colours2 = ['blue','dodgerblue','black','magenta','orange','orange']






fig, ax = plt.subplots(2, 2, sharex='all')#,sharey='all')#, constrained_layout=True)
for m,mode in enumerate(['ldfa/','']):
    for v,vib in enumerate(['v03','v16']):
        time_axis = np.loadtxt(vib+'/'+mode+'timeaxis_1.txt')
        projected_tensors = vjp.read_array(vib+'/'+mode+'projected_tensors1.txt')
        projected_velocities = vjp.read_array(vib+'/'+mode+'projected_velocities_1.txt')
        projected_velocities = np.square(projected_velocities)

        ax2 = ax[v,m].twinx()
        i = 0
        if m ==0 and v == 0:
            ax[v,m].plot(time_axis,projected_velocities[:,i],label=labels1[i],linestyle='-',**line_args,color=colours1[i])
            ax2.plot(time_axis,projected_tensors[:,i,i],label=labels2[i],linestyle='--',**line_args,color=colours2[i])
        else:
            ax[v,m].plot(time_axis,projected_velocities[:,i],linestyle='-',**line_args,color=colours1[i])
            ax2.plot(time_axis,projected_tensors[:,i,i],linestyle='--',**line_args,color=colours2[i])
        ax[v,m].yaxis.label.set_color('red')
        ax[v,m].tick_params(axis='y',which='both', colors='red')
        ax2.yaxis.label.set_color('blue')
        ax2.tick_params(axis='y',which='both', colors='blue')

        
        ax2.set_ylim(bottom=0.0,top=1.2)
        ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax2.yaxis.set_major_locator(MultipleLocator(0.2))
        ax2.set_ylabel(r'$\Lambda_{ij}$ / ps$^{-1}$',color='blue')#,color='white')





for i in range(2):
    for j in range(2):
        #ax[i,j].tick_params(axis='both', which='major')
        ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i,j].set_xlim(left=0.0,right=0.4)
        ax[i,j].set_ylim(bottom=0.0,top=0.004)
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.0005))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.001))

        ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i,j].xaxis.set_major_locator(MultipleLocator(0.1))
        if i ==1:
            ax[i,j].set_xlabel('Time / ps')#,color='white')
        ax[i,j].set_ylabel(r'$\mid\dot{R}\mid ^{2}$ / $\mathrm{ \AA{}}^{2}$ fs$^{-2}$',color='red')#,color='white')


annotate_args = {'xy' : (0.03,0.89), 'xycoords' : 'axes fraction'}
ax[0,0].annotate('(a)',ha="left", **annotate_args)
ax[0,1].annotate('(b)',ha="left", **annotate_args)
ax[1,0].annotate('(c)',ha="left", **annotate_args)
ax[1,1].annotate('(d)',ha="left", **annotate_args)

annotate_args = {'xy' : (0.97,0.89), 'xycoords' : 'axes fraction'}
ax[0,0].annotate(r'$v_i = 3$',ha="right", **annotate_args)
ax[0,1].annotate(r'$v_i = 3$',ha="right", **annotate_args)
ax[1,0].annotate(r'$v_i = 16$',ha="right", **annotate_args)
ax[1,1].annotate(r'$v_i = 16$',ha="right", **annotate_args)

annotate_args = {'xy' : (0.97,0.75), 'xycoords' : 'axes fraction'}
ax[0,1].annotate(r'ODF',ha="right", **annotate_args)
ax[0,0].annotate(r'LDFA',ha="right", **annotate_args)
ax[1,1].annotate(r'ODF',ha="right", **annotate_args)
ax[1,0].annotate(r'LDFA',ha="right", **annotate_args)



plt.subplots_adjust(wspace=0.7)





fig.set_figheight(4.0)
fig.set_figwidth(6.5)

#plt.gcf().subplots_adjust(left=0.3,bottom=0.3,right=0.8)
fig.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.), loc='center')


fig.savefig('figs13.pdf',transparent=True,bbox_inches='tight')
fig.savefig('figs13.png',transparent=True,bbox_inches='tight',dpi=300)
fig.savefig('figs13.eps',transparent=True,bbox_inches='tight')
fig.savefig('figs13.tiff',transparent=True,bbox_inches='tight',dpi=600)