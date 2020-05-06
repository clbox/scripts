import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

line_settings = {
    'marker' : 's',
    'alpha' : 0.6,
    'markersize' : 4
}
fontsize=15
def plot_settings(ax):
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.legend(fontsize=fontsize,fancybox=True,framealpha=0,loc=0)
    ax.set_xlim(0)
    ax.set_ylim(0)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(0.025))
    return

filename = sys.argv[1]


n_files = len(sys.argv[2:])
print('number of files: ' + str(n_files))
test = np.loadtxt(sys.argv[2]+filename,delimiter='  ')
dim = np.shape(test)[0]
time_step = 1 

all_tensors = np.zeros((n_files,dim,dim))
time_axis = np.zeros((n_files))
for i,output_dir in enumerate(sys.argv[2:]):

    all_tensors[i,:,:] = np.loadtxt(output_dir+filename,delimiter='  ')
    time_axis[i] = int(output_dir.replace('/',''))*time_step

#FIG: RELAXATION RATES
fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')
#labels = [r'$d$',r'$\phi$',r'$\theta$',r'$X$',r'$Y$',r'$Z$']
labels = [r'$d$',r'$\theta$',r'$\phi$',r'$X$',r'$Y$',r'$Z$']
for i in range(dim):
    ax.plot(time_axis,all_tensors[:,i,i],label=labels[i],**line_settings)
plot_settings(ax)
fig.set_figheight(15*0.393701)
fig.set_figwidth(15*0.393701)
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=20)
fig.text(-0.05, 0.5, r'Relaxation  rate / $\mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=20)
fig.savefig('relaxation_rates.pdf',transparent=True,bbox_inches='tight')

#FIG: FULL PROJECTED TENSOR
fig, ax = plt.subplots(dim, dim, sharex='all', sharey='all')

for i in range(dim):
    for j in range(i,dim):
        ax[i,j].plot(time_axis,all_tensors[:,i,j])
fig.savefig('tensors_plot.pdf',transparent=True,bbox_inches='tight')