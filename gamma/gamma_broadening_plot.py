import gamma.friction_coupling as fc
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

line_width = 0.6
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.6
matplotlib.rcParams['lines.linewidth'] = 0.6

n_mol = 1
markers = ['o','s','^']*n_mol+['.','>','v']*n_mol
colours = ['red','navy','mediumorchid']*n_mol+['maroon','dodgerblue','gold']*n_mol
linestyles = ['-','-','-']*n_mol+['--','--','--']*n_mol

labels = ['\mathrm{O}_x','\mathrm{O}_y','\mathrm{O}_z']*n_mol+['\mathrm{C}_x','\mathrm{C}_y','\mathrm{C}_z']*n_mol

aims1 = 'aims.out'
gamma_files1 = glob.glob("*gamma*k*out")
gamma_files1.sort()

a = fc.friction_gamma_parser(aims_file=aims1,gamma_files=gamma_files1)

sigmas = np.linspace(0.01,1,100)

tensors = []
for sigma in sigmas:
    b = fc.friction_tensor(a,300,sigma,nspin=1)
    tensor = b.calc_tensor()
    tensors.append(tensor)

    dimension = np.shape(tensor)[0]

fig, ax = plt.subplots(1, 1)

tensors = np.array(tensors)
b=0
for i in range(dimension):

    for j in range(dimension):

        if i==j:


            #if i in [0,1,2,12,13,14]:
            ax.plot(sigmas,tensors[:,i,j],color=colours[b],linestyle=linestyles[b],
                label=r'$\Lambda_{{{}{}}}$'.format(labels[i],labels[j]))
                    #mfc='none',markersize=4,marker=markers[b])
            np.savetxt('gamma_diagonal_'+str(i)+'.txt',tensors[:,i,j])
            #else:
                #ax.plot(sigmas,tensors[:,i,j],color=colours[b],linestyle=linestyles[b])

            b += 1


ax.set_xlim(left=0.0)
ax.set_ylim(bottom=0.0)

ax.set_ylabel(r'$\Lambda_{\mathrm{ij}}$ / ps$^{-1}$',color='black')
ax.set_xlabel(r'$\sigma$ / eV',color='black')

ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.legend(fancybox=True,framealpha=0,ncol=2)

fig.set_figheight(4)
fig.set_figwidth(5)
fig.savefig('gamma_broadening.pdf',transparent=True,bbox_inches='tight')
