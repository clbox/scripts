import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)

line_width = 0.8
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.8
matplotlib.rcParams['lines.linewidth'] = 0.8



element = 2

filename_suffix = '/gamma_diagonal_'+str(element)+'.txt'

filename_prefixes = ['6','8','12'] #,'20','32']

sigmas = np.linspace(0.01,1,100)

fig, ax = plt.subplots(1, 1)


for i,k_grid in enumerate(filename_prefixes):
    print(k_grid)
    n_k = int(k_grid)

    y_data = np.loadtxt(k_grid+filename_suffix)

    ax.plot(sigmas,y_data,label=r'${}\times{}\times1$'.format(n_k,n_k))




ax.set_ylabel(r'$\Lambda_{ij}$ / ps$^{-1}$',color='black')
ax.set_xlabel(r'$\sigma$ / eV',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.legend(fancybox=True,framealpha=0,ncol=2)

ax.set_xlim(0,1)
ax.set_ylim(0,0.7)

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))

ax.text(0.99,0.5,r'$2\times$ c($2\times2$) CO/Cu(100)', ha='right',fontsize=12)
ax.text(0.99,0.46,r'6 Cu layers', ha='right',fontsize=12)


fig.set_figheight(4)
fig.set_figwidth(5)
fig.savefig('gamma_broadenings'+str(element)+'.pdf',transparent=True,bbox_inches='tight')
