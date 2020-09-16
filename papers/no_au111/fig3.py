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
pes_args = {'marker' : 'v', 'linestyle' : 'None','color' : 'green', 'label' : r'ODF PES(2)', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'EXPT', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : '#F5C799', 'markerfacecolor' : 'white', 'label' : r'MDEF Ref', 'alpha' : 1.0}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : '#9ABD8F', 'markerfacecolor' : 'white', 'label' : r'IESH Ref', 'alpha' : 1.0}
d4_args = {'marker' : '^', 'linestyle' : '--','color' : 'indigo', 'label' : r'ODF [$\Lambda_{rr} \times 4$]', 'alpha' : 1.0}
annotate_args = {'xy' : (0.03,0.91), 'xycoords' : 'axes fraction'}

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
v3to3_err = np.array([0.829,0.737,0.875,0.612,0.500,0.493,0.336,0.454])-exp_v3tov3[:,1]
v3to2_err = np.array([0.151,0.270,0.309,0.395,0.454,0.507,0.520,0.539])-exp_v3tov2[:,1]
v3to1_err = np.array([0.0987,0.0724,0.0592,0.118,0.184,0.191,0.237,0.243])-exp_v3tov1[:,1]
errs = [v3to1_err,v3to2_err,v3to3_err]
v3_exp = [exp_v3tov1,exp_v3tov2,exp_v3tov3]


#All bounce
# iesh3to3_pos = np.array([0.356,0.401,0.461,0.531,0.572,0.591,0.632,0.638,0.616,0.677,0.670,0.673])
# iesh3to2_pos = np.array([0.294,0.345,0.371,0.336,0.317,0.342,0.329,0.332,0.358,0.317,0.320,0.307])
# iesh3to1_pos = np.array([0.483,0.379,0.274,0.236,0.207,0.156,0.122,0.112,0.109,0.0835,0.0899,0.103])
#iesh_pos = [iesh3to1_pos,iesh3to2_pos,iesh3to3_pos]

# ef3to3_pos = np.array([0.369,0.207,0.166,0.204,0.229,0.277,0.293,0.305,0.353,0.347,0.369,0.397])
# ef3to2_pos = np.array([0.421,0.615,0.650,0.650,0.685,0.694,0.707,0.726,0.697,0.707,0.688,0.662])
# ef3to1_pos = np.array([0.407,0.309,0.299,0.252,0.182,0.122,0.0835,0.0359,0.0137,0.0105,0.00736,0.00419])
# ef_pos = [ef3to1_pos,ef3to2_pos,ef3to3_pos]




fig, ax = plt.subplots(2, 2)#, sharex='all',sharey='all')#, constrained_layout=True)

final_state = 1
for i in range(2):
    for j in range(2):
        if i == 0 and j == 0:
            exp_v2tov1 = np.loadtxt('v02_exp_1.txt',delimiter=',')
            a = ax[0,0].plot(exp_v2tov1[:,0],exp_v2tov1[:,1],markersize=4,**exp_args)
        else:
            directory =  './v03/translational/'
            ef_results = np.loadtxt(directory+'ef_'+str(final_state)+'_SB.txt',delimiter=',')
            iesh_results = np.loadtxt(directory+'iesh_'+str(final_state)+'_SB.txt',delimiter=',')
            iesh_pos = np.loadtxt(directory+'iesh_'+str(final_state)+'_SB_pos.txt',delimiter=',')
            ef_pos = np.loadtxt(directory+'ef_'+str(final_state)+'_SB_pos.txt',delimiter=',')

            err = iesh_pos - iesh_results[:,1]
            a = ax[i,j].errorbar(iesh_results[:,0],iesh_results[:,1],markersize=4,
                yerr=err,capsize=3,elinewidth=1,zorder=-11,**iesh_args)

            err = ef_pos - ef_results[:,1] 
            a = ax[i,j].errorbar(ef_results[:,0],ef_results[:,1],markersize=4,
                    yerr=err,capsize=3,elinewidth=1,zorder=-10,**ef_args)


            exp = v3_exp[final_state-1]
            err = errs[final_state-1]


            a = ax[i,j].errorbar(exp[:,0],exp[:,1],markersize=4,
                yerr=err,capsize=3,elinewidth=1,zorder=-9,**exp_args)
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

        if 'd4' in os.path.abspath(filename):
            mode = 'd4'

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


map_plot = [(0,0),(0,1),(1,0),(1,1)]
c=0
for initial_state in [2,3]:
    for final_state in [1,2,3]:
        if initial_state == 2 and final_state >= 2:
            continue
        
        for mode in ['ldfa','bomd','tdpt']: #,'pes','d4']:
            if mode=='tdpt':
                mode_args = tdpt_args.copy()
                zorder=3
            if mode=='bomd':
                mode_args = bomd_args.copy()
                zorder=1
            if mode=='ldfa':
                mode_args = ldfa_args.copy()
                zorder=2
            if mode=='d4':
                mode_args = d4_args.copy()

            idx = np.argwhere((all_modes==mode) & (all_initials == initial_state) & (all_finals == final_state)).flatten()
            incidence_es = all_eis[idx]
            ratios = all_ratios[idx]

            order = np.argsort(incidence_es)
            incidence_es = incidence_es[order]
            ratios = ratios[order]
            a = ax[map_plot[c]].plot(incidence_es,ratios,**mode_args,markersize=4,markeredgecolor='black',zorder=zorder)
        ax[map_plot[c]].annotate(r'$\nu_i=$'+str(initial_state)+r'$\rightarrow$'+r'$\nu_f=$'+str(final_state),xy=(0.2,0.91),xycoords='axes fraction')
        #ax[map_plot[c]].annotate(str(initial_state)+r'$\rightarrow$'+str(final_state),xy=(0.5,0.9),xycoords='axes fraction')

        c+=1

    

# 
###########################


letters = [r'(a)',r'(b)',r'(c)',r'(d)']
c=0
for i in range(2):
    for j in range(2):
        ax[i,j].annotate(letters[c],ha="left", **annotate_args)
        c+=1
        if i == 0 and j ==0:
            ax[i,j].tick_params(axis='both', which='major')
            ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
            ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.05))
            ax[i,j].xaxis.set_major_locator(MultipleLocator(0.3))
            continue
        font='Arial'
        for tick in ax[i,j].get_xticklabels():
            tick.set_fontname(font)
        for tick in ax[i,j].get_yticklabels():
            tick.set_fontname(font)

        ax[i,j].set_xlim(0,1.2)
        ax[i,j].tick_params(axis='both', which='major')
        ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i,j].xaxis.set_major_locator(MultipleLocator(0.3))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.05))

        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.2))

ax[0,1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[0,0].set_xlabel(r'$E_i$ / eV',fontname=font,color='black')
ax[0,1].set_xlabel(r'$E_i$ / eV',fontname=font,color='black')
ax[1,0].set_xlabel(r'$E_i$ / eV',fontname=font,color='black')
ax[1,1].set_xlabel(r'$E_i$ / eV',fontname=font,color='black')

ax[0,0].set_ylabel(r'$P(1)\ /\ P(2)$',fontname=font,color='black')
# ax[0,1].set_ylabel(r'$B(\nu_f)$',fontname=font,color='black')
# ax[1,0].set_ylabel(r'$B(\nu_f)$',fontname=font,color='black')

ax[0,1].set_ylabel(r'$P(\nu_f)/[P(1)+P(2)+P(3)]$',fontname=font,color='black')
ax[1,0].set_ylabel(r'$P(\nu_f)/[P(1)+P(2)+P(3)]$',fontname=font,color='black')
#ax[0,1].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
ax[1,1].yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

ax[0,1].set_ylim(0,0.3)
ax[1,0].set_ylim(0,1.2)
ax[1,1].set_ylim(0,1.2)

ax[0,0].set_xlim(0,0.7)
ax[0,0].set_yscale('log')  
labels= (1e-4,1e-3,1e-2,1e-1,1)
ax[0,0].set_yticks(labels)
#ax[0,0].set_yticklabels(labels)
ax[0,0].set_ylim(1e-3,1)


ax[0,1].yaxis.set_label_position("right")
ax[0,1].yaxis.tick_right()
# ax[1,1].yaxis.set_label_position("right")
# ax[1,1].yaxis.tick_right()

labels= (0,0.2,0.4,0.6,0.8,1.0)
ax[1,0].set_yticks(labels)

handles,labels = ax[1,1].get_legend_handles_labels()

handles = [handles[1], handles[3], 
            handles[0], handles[4], 
            handles[2], handles[5]]
labels = [labels[1], labels[3], 
            labels[0],labels[4], 
            labels[2], labels[5]]

fig.set_figheight(4.)
fig.set_figwidth(3.25)
plt.legend(handles=handles,labels=labels,ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(-0.0, 2.5), loc='center')
plt.subplots_adjust(hspace=0.3,wspace=0.2)
fig.savefig('fig3.pdf',transparent=True,bbox_inches='tight',dpi=300)
fig.savefig('fig3.eps',transparent=False,bbox_inches='tight')
