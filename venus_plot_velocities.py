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

time_axis = np.loadtxt(sys.argv[1])

projected_velocities = vjp.read_array(sys.argv[2])

labels = [r'$\Lambda_{rr}$',r'$\Lambda_{\theta\theta}$',r'$\Lambda_{\phi\phi}$',r'$\Lambda_{XX}$',r'$\Lambda_{YY}$',r'$\Lambda_{ZZ}$']
line_args = {'linewidth' : 1.0}
colours = ['red','dodgerblue','black','magenta','orange','orange']


ndim = np.shape(projected_velocities)[1]


projected_velocities = np.square(projected_velocities)

fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)
for i in range(ndim):
    if i in [3,4,5]:
        continue
    if i <=2:
        linestyle = '-'
    else:
        linestyle='-'
    ax.plot(time_axis,projected_velocities[:,i],label=labels[i],linestyle=linestyle,**line_args,color=colours[i])


font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)

ax.tick_params(axis='both', which='major', labelsize=12)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=15)

ax.yaxis.set_minor_locator(MultipleLocator(0.0005))
ax.yaxis.set_major_locator(MultipleLocator(0.001))

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(0.1))

ax.set_xlim(left=0.0,right=0.4)
ax.set_ylim(bottom=0.0,top=0.005)

fig.set_figheight(2.0)
fig.set_figwidth(3.25)

plt.gcf().subplots_adjust(left=0.3,bottom=0.3)


ax.set_xlabel('Time / ps',fontsize=12,fontname=font)#,color='white')
ax.set_ylabel(r'$\mid\dot{R}\mid ^{2}$ / $\mathrm{ \AA{}}^{2}$ fs$^{-2}$',fontsize=12,fontname=font)#,color='white')
fig.savefig('velocities.pdf',transparent=True)#,bbox_inches='tight')


fig.legend(ncol=6,fontsize=12,fancybox=True,framealpha=0)
fig.set_figwidth(7)
ax.remove()
fig.savefig('legend.pdf',transparent=True)#,bbox_inches='tight')
