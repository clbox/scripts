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
import itertools
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
bomd_pes_args = {'marker' : 'v', 'linestyle' : '-','color' : 'cyan', 'label' : r'BOMD[RS]', 'alpha' : 1.0}
odf_pes_args = {'marker' : 'D', 'linestyle' : '-.','color' : 'orange', 'label' : r'ODF[RS]', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'Expt', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : '#F5C799', 'markerfacecolor' : 'white', 'label' : r'MDEF Ref', 'alpha' : 1.0}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : '#9ABD8F', 'label' : r'IESH Ref', 'alpha' : 1.0}
mfcs = ['None','']
annotate_args = {'xy' : (0.03,0.91), 'xycoords' : 'axes fraction'}

results = {'initial' : [], 'final' : [], 'mode' : [], 'incidence_e' : [], 'ratio' : []}




fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)


exp1 = np.loadtxt('exp/1/1050.txt',delimiter=',')
exp2 = np.loadtxt('exp/2/1050.txt',delimiter=',')
a = ax.plot(exp1[:,0],exp1[:,1,],color='black',marker='o',linestyle='-',markersize=4,markeredgecolor='black',label='Expt')
a = ax.plot(exp2[:,0],exp2[:,1,],mfc='None',color='black',marker='s',linestyle='--',markersize=4,markeredgecolor='black')
exp1 = np.loadtxt('ref/1/1050.txt',delimiter=',')
exp2 = np.loadtxt('ref/2/1050.txt',delimiter=',')
a = ax.plot(exp1[:,0],exp1[:,1,],mfc='#9ABD8F',markersize=4,markeredgecolor='black',**iesh_args)
a = ax.plot(exp2[:,0],exp2[:,1,],mfc='white',markersize=4,markeredgecolor='black',**iesh_args)

for i,filename in enumerate(filenames):
    ei = float(os.path.basename(os.path.dirname(filename)))
    dis = np.loadtxt(filename)
    mode_args = None

    if len(np.shape(dis))==1:
        dis = dis.reshape(1,2)
    states = dis[:,0]
    pops = dis[:,1]


    if 'tdpt' in os.path.abspath(filename):
        mode = 'tdpt'

    elif 'bomd' in os.path.abspath(filename):
        mode = 'bomd'

    elif 'ldfa' in os.path.abspath(filename):
        mode = 'ldfa'

    if 'pes' in os.path.abspath(filename):
        mode = 'pes'
        if 'tdpt' in os.path.abspath(filename):
            mode = 'odf[rs]'
        if 'bomd' in os.path.abspath(filename):
            mode = 'bomd[rs]'



    if 'v00' in os.path.abspath(filename):
        initial=2
        final_states = [1,2]
        
        idx = np.argwhere(states==0)
        total = pops[idx]


    for f,state in enumerate(states):
        if state in final_states:
            final = state
            ratio = pops[f]

            results['initial'].append(initial)
            results['final'].append(final)
            results['mode'].append(mode)
            results['incidence_e'].append(ei)
            results['ratio'].append(ratio)

all_modes = np.array(results['mode'])
all_eis = np.array(results['incidence_e'])
all_ratios = np.array(results['ratio'])
all_initials = np.array(results['initial'])
all_finals = np.array(results['final'])

initial_state = 0
for f,final_state in enumerate(final_states):
    for mode in ['bomd','tdpt','bomd[rs]','odf[rs]']:
        if mode=='tdpt':
            mode_args = tdpt_args.copy()
            zorder=3
            model = 'ODF'
        if mode=='bomd':
            mode_args = bomd_args.copy()
            zorder=1
            model = 'BOMD'
        if mode=='ldfa':
            mode_args = ldfa_args.copy()
            zorder=2
            model = 'LDFA'
        if mode=='d4':
            mode_args = d4_args.copy()
            model = 'ODF(rr)*4'
        if mode=='odf[rs]':
            mode_args = odf_pes_args
            model = 'ODF[RS]'
        if mode=='bomd[rs]':
            mode_args = bomd_pes_args
            model = 'BOMD[RS]'

        idx = np.argwhere((all_modes==mode) & (all_finals == final_state)).flatten()
        incidence_es = all_eis[idx]
        ratios = all_ratios[idx]

        order = np.argsort(incidence_es)
        incidence_es = incidence_es[order]
        ratios = ratios[order]
        if f == 0:
            a = ax.plot(incidence_es,ratios,**mode_args,markersize=4,markeredgecolor='black',zorder=zorder)
        elif f == 1:
            a = ax.plot(incidence_es,ratios,**mode_args,mfc='None',markersize=4,markeredgecolor='black',zorder=zorder)

ax.annotate(r'$v_i=$'+str(initial_state)+r'$\rightarrow$'+r'$v_f=$'+str(1),xy=(0.3,0.91),xycoords='axes fraction')
ax.annotate(r'$v_i=$'+str(initial_state)+r'$\rightarrow$'+r'$v_f=$'+str(2),xy=(0.2,0.36),xycoords='axes fraction')
        #ax[map_plot[c]].annotate(str(initial_state)+r'$\rightarrow$'+str(final_state),xy=(0.5,0.9),xycoords='axes fraction')
    
ax.set_yscale('log')  
ax.set_ylim(10e-6,1e-1)
ax.set_xlim(0,1200)
###########################

ax.set_xlabel('$T_s$ / K')
ax.set_ylabel('$P(v_f)$')

handles,labels = ax.get_legend_handles_labels()
print(labels)
handles = handles[0:2] + handles[3:7]
labels = labels[0:2] + labels[3:7]

fig.set_figheight(3.)
fig.set_figwidth(3.25)
plt.legend(handles=handles,labels=labels,ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 1.1), loc='center')
plt.subplots_adjust(hspace=0.3,wspace=0.2)
fig.savefig('excitation.pdf',transparent=True,bbox_inches='tight')
