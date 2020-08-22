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
annotate=True
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"

filenames = sys.argv[1:]



tdpt_args = {'marker' : 'o', 'linestyle' : '--','color' : 'mediumorchid', 'label' : r'ODF', 'alpha' : 1.0}
d4_args = {'marker' : '^', 'linestyle' : '--','color' : 'indigo', 'label' : r'ODF (r) $\times 4$', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'EXPT', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : 'darkorange', 'markerfacecolor' : 'white', 'label' : r'EF ref', 'alpha' : 0.5}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : 'green', 'markerfacecolor' : 'white', 'label' : r'IESH ref', 'alpha' : 0.5}
tdpt_pes_args = {'marker' : '>', 'linestyle' : '-','color' : 'orange', 'label' : r'ODF PES$_\mathrm{rs}$', 'alpha' : 1.0}
bomd_pes_args = {'marker' : 'v', 'linestyle' : '-','color' : 'green', 'label' : r'BOMD PES$_\mathrm{rs}$', 'alpha' : 1.0}
annotate_args = {'xy' : (0.98,0.9), 'xycoords' : 'axes fraction'}

results = {'mode' : [], 'incidence_es' : [], 'ratios' : []}

fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)


exp = np.loadtxt('v02_trapped.txt') #curve
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(func, exp[:,0], exp[:,1])

exp = np.loadtxt('/Users/u1865573/work/colab/no_au111/venus/analysis/v02/translational/exp_trapped.txt')

#reads in multiple absolute_pop.txt 
for i,filename in enumerate(filenames):
    ei = float(os.path.basename(os.path.dirname(filename)))
    dis = np.loadtxt(filename)
    mode_args = None

    

    if -1 not in dis[:,0]:
        continue

    ratio = dis[-1,1]/np.sum(dis[:,1])
    
    if 'tdpt' in os.path.abspath(filename):
        mode = 'tdpt'

    if 'bomd' in os.path.abspath(filename):
        mode = 'bomd'

    if 'ldfa' in os.path.abspath(filename):
        mode = 'ldfa'

    if 'd4' in os.path.abspath(filename):
        mode = 'd4'

    if 'pes' in os.path.abspath(filename):
        if 'tdpt' in os.path.abspath(filename):
            mode = 'tdpt_pes'

        if 'bomd' in os.path.abspath(filename):
            mode = 'bomd_pes'



    if ei not in[97,300,420,425,640]:
        continue
        
    # if plot_v02:
    #     if ei not in[97,300,420,425,640]:
    #         continue
            
    #     if mode == 'bomd' and ei <= 100:
    #         continue

    results['mode'].append(mode)
    results['incidence_es'].append(ei/1000)
    results['ratios'].append(ratio)

all_modes = np.array(results['mode'])
all_eis = np.array(results['incidence_es'])
all_ratios = np.array(results['ratios'])

for mode in ['bomd','ldfa','tdpt','d4','tdpt_pes','bomd_pes']:
    print(mode)
    idx = np.argwhere(all_modes==mode)

    idx = idx[:,0]


    
    incidence_es = all_eis[idx]
    print(incidence_es)
    ratios = all_ratios[idx]

    order = np.argsort(incidence_es)
    incidence_es = incidence_es[order]
    ratios = ratios[order]

    if mode=='tdpt':
        mode_args = tdpt_args.copy()
    if mode=='bomd':
        mode_args = bomd_args.copy()
    if mode=='ldfa':
        mode_args = ldfa_args.copy()
    if mode=='d4':
        mode_args = d4_args.copy()
    if mode=='tdpt_pes':
        mode_args = tdpt_pes_args.copy()
    if mode=='bomd_pes':
        mode_args = bomd_pes_args.copy()

    mode_args['linestyle'] = 'None'
    a = ax.plot(incidence_es,ratios,**mode_args,markersize=6,markeredgecolor='black')


ax.plot(np.linspace(0,1,100), func(np.linspace(0,1,100), *popt), '--', color='black')
ax.errorbar(exp[:,0],exp[:,1],yerr=exp[:,2]-exp[:,1],markersize=4,capsize=3,elinewidth=1,zorder=-10,linestyle='none',color='black')
#Fake line for legend entry

ax.errorbar([-100,-50],[-60,-70],yerr=10,markersize=4,capsize=3,elinewidth=1,zorder=-10,linestyle='--',color='black',label='EXPT')
###########################
font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)

ax.annotate(r'$\nu_i = 2$',ha="right", **annotate_args)

ax.tick_params(axis='both', which='major')
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=15)

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(0.2))

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.2))

ax.set_xlim(0,0.8)
ax.set_ylim(0,1)


#ax.set_yscale('log')  
fig.set_figheight(3.0)
fig.set_figwidth(3.25)


#fig.set_constrained_layout_pads(w_pad=0, h_pad=0)

ax.set_xlabel('Incidence energy / eV',fontname=font)#,color='white')
ax.set_ylabel('Population',fontname=font)#,color='white')

#plt.gcf().subplots_adjust(left=0.3,bottom=0.3)

plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.1), loc='center')
plt.subplots_adjust(hspace=0.3,wspace=0.2)
fig.savefig('fig6.pdf',transparent=True,bbox_inches='tight')
fig.savefig('fig6.eps',transparent=False,bbox_inches='tight')
