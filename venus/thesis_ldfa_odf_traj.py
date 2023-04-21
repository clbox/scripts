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
plt.style.use('clb_publication_thesis')

time_odf_file = "/Users/u1865573/work/colab/no_au111/venus/analysis/v03/overview/tdpt/eft/plot_test/v03/timeaxis_1.txt"
tens_odf_file = "/Users/u1865573/work/colab/no_au111/venus/analysis/v03/overview/tdpt/eft/plot_test/v03/projected_tensors1.txt"
vel_odf_file = "/Users/u1865573/work/colab/no_au111/venus/analysis/v03/overview/tdpt/eft/plot_test/v03/projected_velocities_1.txt"


time_ldfa_file = "/Users/u1865573/work/colab/no_au111/venus/analysis/v03/overview/tdpt/eft/plot_test/v03/ldfa/timeaxis_1.txt"
tens_ldfa_file = "/Users/u1865573/work/colab/no_au111/venus/analysis/v03/overview/tdpt/eft/plot_test/v03/ldfa/projected_tensors1.txt"
vel_ldfa_file = "/Users/u1865573/work/colab/no_au111/venus/analysis/v03/overview/tdpt/eft/plot_test/v03/ldfa/projected_velocities_1.txt"

time_axis_odf = np.loadtxt(time_odf_file)
projected_tensors_odf = vjp.read_array(tens_odf_file)
projected_velocities_odf = vjp.read_array(vel_odf_file)


time_axis_ldfa = np.loadtxt(time_ldfa_file)
projected_tensors_ldfa = vjp.read_array(tens_ldfa_file)
projected_velocities_ldfa = vjp.read_array(vel_ldfa_file)

projected_velocities_odf = np.square(projected_velocities_odf)
projected_velocities_ldfa = np.square(projected_velocities_ldfa)

time_axes = [time_axis_ldfa,time_axis_odf]
proj_tensors = [projected_tensors_ldfa,projected_tensors_odf]
proj_vels = [projected_velocities_ldfa,projected_velocities_odf]

labels2 = [r'$\Lambda_{rr}$',r'$\Lambda_{\theta\theta}$',r'$\Lambda_{\phi\phi}$',r'$\Lambda_{XX}$',r'$\Lambda_{YY}$',r'$\Lambda_{ZZ}$']
labels1 = [r'$\dot{R}_{rr}$',r'$\dot{R}_{\theta\theta}$',r'$\dot{R}_{\phi\phi}$',r'$\dot{R}_{XX}$',r'$\dot{R}{YY}$',r'$\dot{R}_{ZZ}$']



ndim = np.shape(projected_velocities_odf)[1]

fig, ax1 = plt.subplots(2, 1, sharex='all',sharey='all')#, constrained_layout=True)



for m in range(2):
    ax2 = ax1[m].twinx()
    for i in range(ndim):
        if i in [1,2,3,4,5]:
            continue
        if i <=2:
            linestyle = '-'
        else:
            linestyle='-'
        ax1[m].plot(time_axes[m],proj_vels[m][:,i],label=labels1[i],marker='None',linewidth=.5,linestyle='-',color='#1F90FF')
        ax2.plot(time_axes[m],proj_tensors[m][:,i,i],label=labels2[i],marker='None',linewidth=1,linestyle='-')

    for ax in [ax1[m],ax2]:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))


    ax1[m].yaxis.set_minor_locator(MultipleLocator(0.0005))
    ax1[m].yaxis.set_major_locator(MultipleLocator(0.001))

    ax1[m].xaxis.set_minor_locator(MultipleLocator(0.05))
    ax1[m].xaxis.set_major_locator(MultipleLocator(0.1))

    ax1[m].set_xlim(left=0.0,right=0.4)
    ax1[m].set_ylim(bottom=0.0,top=0.001)


    ax2.set_ylim(bottom=0.0,top=0.8)
    ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))

    ax2.set_ylabel(r'$\Lambda_{rr}$ / ps$^{-1}$',color='#FF4500')#,color='white')

fig.set_figheight(3)
fig.set_figwidth(1.5748)

# plt.gcf().subplots_adjust(left=0.3,bottom=0.3,right=0.8)


ax1[1].set_xlabel('Time / ps')#,color='white')
ax1[0].set_ylabel(r'$\mid\dot{R_{rr}}\mid ^{2}$ / $\mathrm{ \AA{}}^{2}$ fs$^{-2}$',color='#1F90FF')#,color='white')
ax1[1].set_ylabel(r'$\mid\dot{R_{rr}}\mid ^{2}$ / $\mathrm{ \AA{}}^{2}$ fs$^{-2}$',color='#1F90FF')#,color='white')

fig.savefig('both_new.pdf',transparent=True,bbox_inches='tight')


fig.legend(ncol=6)
fig.set_figwidth(7)
ax1[0].remove()
ax1[1].remove()
ax2.remove()
fig.savefig('legend_new.pdf',transparent=True)#,bbox_inches='tight')