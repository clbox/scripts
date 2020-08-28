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
from scipy.optimize import curve_fit
SMALL_SIZE = 9.5
MEDIUM_SIZE = 9.5
BIGGER_SIZE = 9.5

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
annotate=True
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"

filenames = sys.argv[1:]

colours = ['indigo','maroon','darkgreen','darkorange','navy']

tdpt_args = {'marker' : 'o', 'linestyle' : '-','color' : 'purple', 'label' : r'ODF', 'alpha' : 1.0}
d4_args = {'marker' : '^', 'linestyle' : '--','color' : 'mediumorchid', 'label' : r'ODF (r) $\times 4$', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'Exp', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : 'darkorange', 'markerfacecolor' : 'white', 'label' : r'EF ref', 'alpha' : 0.5}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : 'green', 'markerfacecolor' : 'white', 'label' : r'IESH ref', 'alpha' : 0.5}

annotate_args = {'xy' : (0.72,0.8), 'xycoords' : 'figure fraction'}

results = {'mode' : [], 'incidence_es' : [], 'ratios' : []}

fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)

#reads in multiple absolute_pop.txt 
for i,filename in enumerate(filenames):
    dis = np.loadtxt(filename)
    mode_args = None

    

    if -1 not in dis[:,0]:
        continue

    ratio = dis[-1,1]/np.sum(dis[:,1])
    
    if 'tdpt' in os.path.abspath(filename):
        mode = 'ODF'

    if 'bomd' in os.path.abspath(filename):
        mode = 'BOMD'

    if 'ldfa' in os.path.abspath(filename):
        mode = 'LDFA'

    if 'd4' in os.path.abspath(filename):
        mode = r'ODF [$\Lambda_{rr} \times 4$]'
    
    if 'i4' in os.path.abspath(filename):
        mode += r'[$ \mathbf{\Lambda} \times 4$]'
    if 'pes' in os.path.abspath(filename):
        mode += r'[RS]'

    results['mode'].append(mode)
    results['ratios'].append(ratio)

all_modes = np.array(results['mode'])
all_ratios = np.array(results['ratios'])

# for mode in ['bomd','ldfa','tdpt','d4']:
#     print(mode)
#     idx = np.argwhere(all_modes==mode)

#     idx = idx[:,0]


    
#     incidence_es = all_eis[idx]
#     print(incidence_es)
#     ratios = all_ratios[idx]

#     order = np.argsort(incidence_es)
#     incidence_es = incidence_es[order]
#     ratios = ratios[order]

#     if mode=='tdpt':
#         mode_args = tdpt_args.copy()
#     if mode=='bomd':
#         mode_args = bomd_args.copy()
#     if mode=='ldfa':
#         mode_args = ldfa_args.copy()
#     if mode=='d4':
#         mode_args = d4_args.copy()

#     mode_args['linestyle'] = 'None'
#     a = ax.plot(incidence_es,ratios,**mode_args,markersize=6,markeredgecolor='black')
print(all_modes)
print(all_ratios)
for i,v in enumerate(['02','03','11','15','16']):
    if v in filenames[0]:
        bar_colour = colours[i]
        vib_state = int(v)
ax.barh(all_modes,all_ratios,color=bar_colour,edgecolor='black')
ax.annotate(r'$\nu_i = {}$'.format(vib_state), **annotate_args)




ax.set_yticklabels(labels=all_modes)
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=15)

# ax.xaxis.set_minor_locator(MultipleLocator(0.05))
# ax.xaxis.set_major_locator(MultipleLocator(0.2))

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(0.1))

#ax.set_xlim(0,0.8)
ax.set_xlim(0,0.5)


#ax.set_yscale('log')  
fig.set_figheight(2.0)
fig.set_figwidth(3.25)


#fig.set_constrained_layout_pads(w_pad=0, h_pad=0)

ax.set_ylabel('Model')#,color='white')
ax.set_xlabel('Population')#,color='white')

plt.gcf().subplots_adjust(left=0.4,bottom=0.2)

#fig.text(0.5, 0.00, r"Final vibrational state ($\nu_f$)", ha='center',fontsize=15)
#fig.text(0.01, 0.5, 'Population', va='center', rotation='vertical',fontsize=15)
fig.savefig('trapped.pdf',transparent=True)#,bbox_inches='tight')
