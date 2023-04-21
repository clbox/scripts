import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

line_width = 0.8
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.5
matplotlib.rcParams['lines.linewidth'] = 1.2


colours = ['orangered','dodgerblue','seagreen','mediumorchid','maroon','dodgerblue','gold']
linestyles = ['--','-',(0, (5, 10)), '-.', ':']
element = 2

filename_suffix = '/06/gamma_diagonal_'+str(element)+'_proper.txt'
#filename_suffix = '/16/gamma_diagonal_'+str(element)+'.txt'

filename_prefixes = ['12','24','48','64']
#filename_prefixes = ['24','48','64']

sigmas = np.linspace(0.01,1,100)

fig, ax = plt.subplots(1, 1)


for i,k_grid in enumerate(filename_prefixes):
    print(k_grid)
    n_k = int(k_grid)

    y_data = np.loadtxt(k_grid+filename_suffix)


    ax.plot(sigmas,y_data,label=r'$N_k = {}$'.format(n_k),color=colours[i],linestyle=linestyles[i])




ax.set_xlabel(r'$\sigma$ / $\mathrm{eV}$',color='black')
#ax.set_ylabel(r'$\Lambda_{\mathrm{C}_z,\mathrm{C}_z}$ / meV ps$^{-1}$ $\mathrm{\AA{}}^{-2}$',color='black')
ax.set_ylabel(r'$\Lambda_{\mathrm{C}_z,\mathrm{C}_z}$ / $\mathrm{meV ps}$ $\mathrm{\AA{}}^{-2}$',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.xaxis.set_major_locator(MultipleLocator(0.2))
ax.xaxis.set_minor_locator(MultipleLocator(0.02))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.yaxis.set_major_locator(MultipleLocator(0.1))

ax.set_xlim(0,0.8)
ax.set_ylim(0.,0.7)

ax.legend(fancybox=True,framealpha=1,edgecolor='black',handletextpad=0.05,borderpad=0.3,handlelength=1.2)

#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.15), loc='center')
fig.set_figheight(3)
fig.set_figwidth(2)

fig.savefig('gamma_broadenings'+str(element)+'.pdf',transparent=True,bbox_inches='tight')
