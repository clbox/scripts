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
from matplotlib import gridspec
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
d4_args = {'marker' : '^', 'linestyle' : '--','color' : 'indigo', 'label' : r'ODF [$\Lambda_{\mathrm{rr}} \times 4$]', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'EXPT', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : 'darkorange', 'markerfacecolor' : 'white', 'label' : r'EF ref', 'alpha' : 0.5}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : 'green', 'markerfacecolor' : 'white', 'label' : r'IESH ref', 'alpha' : 0.5}
tdpt_pes_args = {'marker' : '>', 'linestyle' : '-','color' : 'orange', 'label' : r'ODF[RS]', 'alpha' : 1.0}
bomd_pes_args = {'marker' : 'v', 'linestyle' : '-','color' : 'green', 'label' : r'BOMD[RS]', 'alpha' : 1.0}
annotate_args = {'xy' : (0.98,0.94), 'xycoords' : 'axes fraction'}

results = {'mode' : [], 'incidence_es' : [], 'ratios' : [], 'vi' : []}

fig = plt.figure()
gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1, 1], width_ratios = [0.7,0.3] )

ax0= fig.add_subplot(gs[:, 0])
ax1= fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[1, 1])

ax = [ax0,ax1,ax2]


#fig, ax = plt.subplots(1, 3)#, constrained_layout=True)


exp = np.loadtxt('v02_trapped.txt') #curve
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(func, exp[:,0], exp[:,1])

exp = np.loadtxt('/Users/u1865573/work/colab/no_au111/venus/analysis/v02/translational/exp_trapped.txt')

#reads in multiple absolute_pop.txt 
for i,filename in enumerate(filenames):

    if not 'v02' in filename:
        continue

    v=0
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

    if 'v02' in os.path.abspath(filename):
        v=2

    if 'v11' in os.path.abspath(filename):
        v=11

    if 'v16' in os.path.abspath(filename):
        v=11


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
    results['vi'].append(v)

all_modes = np.array(results['mode'])
all_eis = np.array(results['incidence_es'])
all_ratios = np.array(results['ratios'])
all_v = np.array(results['vi'])

#labels = ['BOMD','LDFA','ODF',r'ODF [$\Lambda_{\mathrm{rr}} \times 4$]',r'ODF PES$_{\mathrm{rs}}$',r'BOMD PES$_{\mathrm{rs}}$']

for i,mode in enumerate(['bomd','ldfa','tdpt','d4','bomd_pes','tdpt_pes']):
    zorder=i
    print(mode)
    v=2
    idx = np.argwhere((all_modes==mode) & (all_v == v))

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
        zorder=10
    if mode=='ldfa':
        mode_args = ldfa_args.copy()
    if mode=='d4':
        mode_args = d4_args.copy()
        continue
    if mode=='tdpt_pes':
        mode_args = tdpt_pes_args.copy()
    if mode=='bomd_pes':
        mode_args = bomd_pes_args.copy()

    mode_args['linestyle'] = 'None'
    a = ax[0].plot(incidence_es,ratios,**mode_args,markersize=6,markeredgecolor='black',zorder=zorder)


ax[0].plot(np.linspace(0,1,100), func(np.linspace(0,1,100), *popt), '--', color='black')
ax[0].errorbar(exp[:,0],exp[:,1],yerr=exp[:,2]-exp[:,1],markersize=4,capsize=3,elinewidth=1,zorder=-10,linestyle='none',color='black')
#Fake line for legend entry

ax[0].errorbar([-100,-50],[-60,-70],yerr=10,markersize=4,capsize=3,elinewidth=1,zorder=-10,linestyle='--',color='black',label='EXPT')
###########################
ax[0].annotate(r'$\nu_i = 2$',ha="right", **annotate_args)
ax[0].xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=15)
ax[0].xaxis.set_minor_locator(MultipleLocator(0.05))
ax[0].xaxis.set_major_locator(MultipleLocator(0.2))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[0].yaxis.set_major_locator(MultipleLocator(0.2))
ax[0].set_xlim(0,0.8)
ax[0].set_ylim(0,1)





high_results = {'mode' : [], 'incidence_es' : [], 'ratios' : [], 'vi' : [], 'pes' : []}
for i,filename in enumerate(filenames):

    if 'v02' in filename:
        continue

    dis = np.loadtxt(filename)
    mode_args = None

    if -1 not in dis[:,0]:
        continue

    ratio = dis[-1,1]/np.sum(dis[:,1])

    if 'v11' in os.path.abspath(filename):
        v=11

    if 'v16' in os.path.abspath(filename):
        v=16
    
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
        #mode += r' PES$_\mathrm{rs}$'
        pes = 2
    else:
        pes = 1

    high_results['mode'].append(mode)
    high_results['ratios'].append(ratio)
    high_results['vi'].append(v)
    high_results['pes'].append(pes)

all_modes = np.array(high_results['mode'])
all_ratios = np.array(high_results['ratios'])
all_v = np.array(high_results['vi'])
all_pes = np.array(high_results['pes'])
print(all_v)
colours = np.array((['maroon','navy'],['salmon','dodgerblue']))
for i,v in enumerate(['11','16']):
    for p, pes in enumerate([1,2]):
        bar_colour = colours[p,i]
        vib_state = int(v)

        idx = np.argwhere((all_v == vib_state) & (all_pes == pes))

        idx = idx[:,0]
        ratios = all_ratios[idx]
        modes = all_modes[idx]

        annotate_args['xy'] = (0.33,0.05)
        ax[i+1].barh(modes,ratios,color=bar_colour,edgecolor='black')
    ax[i+1].annotate(r'$\nu_i = {}$'.format(vib_state), **annotate_args)


for i in [1,2]:
    ax[i].set_xlim(0,0.5)
    
    ax[i].set_yticklabels(labels=all_modes,rotation='vertical', va = 'center')
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.1))

annotate_args['xy'] = (0.05,0.94)
ax[0].annotate(r'(a)',ha="left", **annotate_args)

annotate_args['xy'] = (0.97,0.89)
ax[1].annotate(r'(b)',ha="right", **annotate_args)
ax[2].annotate(r'(c)',ha="right", **annotate_args)


ax[2].set_xlabel('Population')
ax[1].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

#ax.set_yscale('log')  
fig.set_figheight(3.4)
fig.set_figwidth(3.25)


#fig.set_constrained_layout_pads(w_pad=0, h_pad=0)

ax[0].set_xlabel('Incidence energy / eV')#,color='white')
ax[0].set_ylabel(r'$p_{\mathrm{trap}}$')#,color='white')

ax[0].legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.6, 1.1), loc='center')
#plt.gcf().subplots_adjust(left=0.3,bottom=0.3)

plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.1), loc='center')
plt.subplots_adjust(hspace=0.1,wspace=0.3)
fig.savefig('fig6.pdf',transparent=True,bbox_inches='tight')
fig.savefig('fig6.eps',transparent=False,bbox_inches='tight')
