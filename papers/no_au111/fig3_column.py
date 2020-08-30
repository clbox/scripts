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

tdpt_args = {'marker' : 'o', 'linestyle' : '-','color' : 'mediumorchid', 'label' : r'ODF', 'alpha' : 1.0}
pes_args = {'marker' : 'v', 'linestyle' : 'None','color' : 'green', 'label' : r'ODF PES(2)', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'Exp', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : 'darkorange', 'markerfacecolor' : 'white', 'label' : r'EF ref', 'alpha' : 0.5}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : 'green', 'markerfacecolor' : 'white', 'label' : r'IESH ref', 'alpha' : 0.5}

annotate_args = {'xy' : (0.05,0.80), 'xycoords' : 'axes fraction'}

results = {'initial' : [], 'final' : [], 'mode' : [], 'incidence_e' : [], 'ratio' : []}
exp_v3tov1 = np.array([[0.1131, 0.1048],
[0.2625, 0.07947],
[0.3841, 0.07015],
[0.5221, 0.1188],
[0.6494, 0.1930],
[0.8163, 0.1613],
[0.9473, 0.2196],
[1.069, 0.2295]])
exp_v3tov2 = np.array([[0.1036, 0.1572],
[0.2553, 0.2679],
[0.3767, 0.2574],
[0.5076, 0.3680],
[0.6391, 0.4016],
[0.8086, 0.4463],
[0.9366, 0.4909],
[1.055, 0.4475]])
exp_v3tov3 = np.array([[0.09586, 0.7539],
[0.2454, 0.7014],
[0.3667, 0.7088],
[0.5031, 0.5429],
[0.6322, 0.4303],
[0.7988, 0.4178],
[0.9347, 0.3186],
[1.045, 0.3725]])

v3_exp = [exp_v3tov1,exp_v3tov2,exp_v3tov3]

fig = plt.figure()
gs0 = gridspec.GridSpec(3,1,height_ratios=[0.3,0.000005,1])

gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0], hspace=0)
ax0 = fig.add_subplot(gs00[0])

gx = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1], hspace=-0.2)

gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[2:], hspace=0.2)
ax1 = fig.add_subplot(gs01[0])
ax2 = fig.add_subplot(gs01[1])
ax3 = fig.add_subplot(gs01[2])

ax = [ax0,ax1,ax2,ax3]


final_state = 1
for i in range(4):
        if i == 0:
            exp_v2tov1 = np.loadtxt('v02_exp_1.txt',delimiter=',')
            a = ax[0].plot(exp_v2tov1[:,0],exp_v2tov1[:,1],markersize=4,**exp_args)
        else:
            ef_results = np.loadtxt('v03_ef_'+str(final_state)+'.txt',delimiter=',')
            iesh_results = np.loadtxt('v03_iesh_'+str(final_state)+'.txt',delimiter=',')
            a = ax[i].plot(ef_results[:,0],ef_results[:,1],markersize=4,**ef_args)
            a = ax[i].plot(iesh_results[:,0],iesh_results[:,1],markersize=4,**iesh_args)

            exp = v3_exp[final_state-1]
            a = ax[i].plot(exp[:,0],exp[:,1],markersize=4,**exp_args)
            final_state +=1


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

    if 'v02' in os.path.abspath(filename):
        initial=2
        final_states = [1]
        
        idx = np.argwhere(states==2)
        total = pops[idx]


    elif 'v03' in os.path.abspath(filename):
        initial=3
        final_states = [1,2,3]

        total = 0
        for f in range(len(states)):
            if states[f] in final_states:
                total += pops[f]

    for f,state in enumerate(states):
        if state in final_states:
            final = state
            ratio = pops[f]/total

            results['initial'].append(initial)
            results['final'].append(final)
            results['mode'].append(mode)
            results['incidence_e'].append(ei/1000)
            results['ratio'].append(ratio)

all_modes = np.array(results['mode'])
all_eis = np.array(results['incidence_e'])
all_ratios = np.array(results['ratio'])
all_initials = np.array(results['initial'])
all_finals = np.array(results['final'])


map_plot = [(0),(1),(2),(3)]
c=0
for initial_state in [2,3]:
    for final_state in [1,2,3]:
        if initial_state == 2 and final_state >= 2:
            continue
        
        print('initial ' + str(initial_state))
        print('final ' + str(final_state))
        for mode in ['ldfa','bomd','tdpt']:
            if mode=='tdpt':
                mode_args = tdpt_args.copy()
            if mode=='bomd':
                mode_args = bomd_args.copy()
            if mode=='ldfa':
                mode_args = ldfa_args.copy()

            idx = np.argwhere((all_modes==mode) & (all_initials == initial_state) & (all_finals == final_state)).flatten()
            incidence_es = all_eis[idx]
            ratios = all_ratios[idx]

            order = np.argsort(incidence_es)
            incidence_es = incidence_es[order]
            ratios = ratios[order]
            print(c)
            a = ax[map_plot[c]].plot(incidence_es,ratios,**mode_args,markersize=4,markeredgecolor='black')
        ax[map_plot[c]].annotate(r'$\nu_i=$'+str(initial_state)+r'$\rightarrow$'+r'$\nu_f=$'+str(final_state),xy=(0.2,0.8),xycoords='axes fraction')
        #ax[map_plot[c]].annotate(str(initial_state)+r'$\rightarrow$'+str(final_state),xy=(0.5,0.9),xycoords='axes fraction')

        c+=1

    

# 
###########################


letters = [r'(a)',r'(b)',r'(c)',r'(d)']
c=0
for i in range(4):
        ax[i].annotate(letters[c],ha="left", **annotate_args)
        c+=1
        if i == 0:
            ax[i].tick_params(axis='both', which='major')
            ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))
            ax[i].xaxis.set_minor_locator(MultipleLocator(0.05))
            ax[i].xaxis.set_major_locator(MultipleLocator(0.3))
            continue
        font='Arial'
        for tick in ax[i].get_xticklabels():
            tick.set_fontname(font)
        for tick in ax[i].get_yticklabels():
            tick.set_fontname(font)

        ax[i].set_xlim(0,1.2)
        ax[i].tick_params(axis='both', which='major')
        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i].xaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i].xaxis.set_major_locator(MultipleLocator(0.2))
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i].yaxis.set_major_locator(MultipleLocator(0.2))

ax[0].set_xlabel(r'$E_i$ / eV',fontname=font,color='black')
ax[3].set_xlabel(r'$E_i$ / eV',fontname=font,color='black')

ax[0].set_ylabel(r'$P(1)/P(2)$',fontname=font,color='black')
ax[1].set_ylabel(r'$B(\nu_f)$',fontname=font,color='black')
ax[2].set_ylabel(r'$B(\nu_f)$',fontname=font,color='black')

ax[1].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
ax[2].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

ax[1].set_ylim(0,0.6)
ax[2].set_ylim(0,1.)
ax[3].set_ylim(0,1.)

ax[0].set_xlim(0,0.7)
ax[0].set_yscale('log')  
labels= (1e-4,1e-3,1e-2,1e-1,1)
ax[0].set_yticks(labels)
#ax[0,0].set_yticklabels(labels)
ax[0].set_ylim(1e-3,1)


fig.set_figheight(5.)
fig.set_figwidth(3.25)
plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 5.55), loc='center')
plt.subplots_adjust(hspace=0.3)
fig.savefig('fig3.pdf',transparent=True,bbox_inches='tight',dpi=300)
fig.savefig('fig3.eps',transparent=False)#,bbox_inches='tight')
