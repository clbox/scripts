import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
line_settings = {
    'marker' : 's'
}

n_files = len(sys.argv[1:])
print('number of files: ' + str(n_files))
test = np.loadtxt(sys.argv[1]+'projected_tensor.out',delimiter='  ')
dim = np.shape(test)[0]
time_step = 2

labels = [r'$d$',r'$\phi$',r'$\theta$',r'$X$',r'$Y$',r'$Z$']
all_tensors = np.zeros((n_files,dim,dim))
time_axis = np.zeros((n_files))
for i,output_dir in enumerate(sys.argv[1:]):

    all_tensors[i,:,:] = np.loadtxt(output_dir+'projected_tensor.out',delimiter='  ')
    time_axis[i] = int(output_dir.replace('/',''))*time_step

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all')
for i in range(dim):
    ax.plot(time_axis,all_tensors[:,i,i],label=labels[i],**line_settings)
ax.legend()
fig.text(0.5, 0.01, "Time / fs", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'Relaxation  rate / $\mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)

fig.savefig('relaxation_rates.pdf',transparent=True,bbox_inches='tight')




fig, ax = plt.subplots(dim, dim, sharex='all', sharey='all')

for i in range(dim):
    for j in range(i,dim):
        ax[i,j].plot(time_axis,all_tensors[:,i,j])
fig.savefig('projected_tensors.pdf',transparent=True,bbox_inches='tight')