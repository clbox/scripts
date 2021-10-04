import gamma.friction_coupling as fc
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

line_width = 0.4
matplotlib.rcParams['axes.linewidth'] = line_width
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['lines.markeredgewidth'] = 0.6
matplotlib.rcParams['lines.linewidth'] = 0.6

markers = ['o','s','^','.','>','v']
colours = ['red','navy','mediumorchid','maroon','dodgerblue','gold']
linestyles = ['-','-','-','--','--','--']


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

            ax.plot(sigmas,tensors[:,i,j],color=colours[b],linestyle=linestyles[b])
                    #mfc='none',markersize=4,marker=markers[b])

            b += 1





ax.xaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=4, width=line_width, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=2, width=line_width, direction='in', right='on')

ax.legend(fancybox=True,framealpha=0)

fig.set_figheight(4)
fig.set_figwidth(5)
fig.savefig('gamma_broadening.pdf',transparent=True,bbox_inches='tight')