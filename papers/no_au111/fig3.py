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

# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
annotate=True
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"


filenames = sys.argv[1:]


fontsize=10

tdpt_args = {'marker' : 'o', 'linestyle' : '-','color' : 'purple', 'label' : r'ODF', 'alpha' : 1.0}
pes_args = {'marker' : 'v', 'linestyle' : 'None','color' : 'green', 'label' : r'ODF PES(2)', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}
exp_args = {'marker' : 's','linestyle' : '-','color' : 'black', 'markerfacecolor' : 'gold', 'label' : r'Exp', 'alpha' : 1.0}
ef_args = {'marker' : 's','linestyle' : '-','color' : 'darkorange', 'markerfacecolor' : 'white', 'label' : r'EF ref', 'alpha' : 0.5}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : 'green', 'markerfacecolor' : 'white', 'label' : r'IESH ref', 'alpha' : 0.5}

annotate_args = {'fontsize' : fontsize, 'xy' : (0.55,0.8), 'xycoords' : 'figure fraction'}

results = {'initial' : [], 'final' : [], 'mode' : [], 'incidence_es' : [], 'ratios' : []}
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

fig, ax = plt.subplots(2, 2, sharex='all',sharey='all')#, constrained_layout=True)

final_state = 1
for i in range(2):
    for j in range(2):
        if i == 0 and j == 0:
            exp_v2tov1 = np.loadtxt('v02_exp_1.txt',delimiter=',')
            a = ax[0,0].plot(exp_v2tov1[:,0],exp_v2tov1[:,1],**exp_args)
        else:
            ef_results = np.loadtxt('ef_'+str(final_state)+'.txt',delimiter=',')
            iesh_results = np.loadtxt('iesh_'+str(final_state)+'.txt',delimiter=',')
            a = ax[i,j].plot(ef_results[:,0],ef_results[:,1],**ef_args)
            a = ax[i,j].plot(iesh_results[:,0],iesh_results[:,1],**iesh_args)

            exp = v3_exp[final_state-1]
            a = ax[i,j].plot(exp[:,0],exp[:,1],**exp_args)
            final_state +=1


for i,filename in enumerate(filenames):
    ei = float(os.path.basename(os.path.dirname(filename)))
    dis = np.loadtxt(filename)
    mode_args = None

    try:
        index = np.argwhere(dis[:,0]==final_state)
    except:
        dis = dis.reshape(1,2)
        index = [[0]]

    try:
        index = index[0]
    except:
        continue
    ratio_final = (dis[index,1])[0]
    if 'tdpt' in os.path.abspath(filename):
        mode = 'tdpt'

    if 'bomd' in os.path.abspath(filename):
        mode = 'bomd'

    if 'ldfa' in os.path.abspath(filename):
        mode = 'ldfa'

    if 'pes' in os.path.abspath(filename):
        mode = 'pes'


    results['mode'].append(mode)
    results['incidence_es'].append(ei/1000)

    if 'v02' in os.path.abspath(filename):
        plot_v02 = True
        try:
            one_index = np.argwhere(dis[:,0]==1)
            ratio_one = (dis[one_index,1])[0]
        except:
            ratio_one = 0

        try:
            two_index = np.argwhere(dis[:,0]==2)
            ratio_two = (dis[two_index,1])[0]
        except:
            ratio_two = 0

        try:
            three_index = np.argwhere(dis[:,0]==3)
            ratio_three = (dis[three_index,1])[0]
        except:
            ratio_three = 0
        
        if final_state == 1:
            ratio = ratio_one / (ratio_two)
        elif final_state == 3:
             ratio = ratio_three / (ratio_two)
        

    if 'v03' in os.path.abspath(filename):
        plot_v03 = True
        try:
            one_index = np.argwhere(dis[:,0]==1)
            ratio_one = (dis[one_index,1])[0]
        except:
            ratio_one = 0

        try:
            two_index = np.argwhere(dis[:,0]==2)
            ratio_two = (dis[two_index,1])[0]
        except:
            ratio_two = 0

        try:
            three_index = np.argwhere(dis[:,0]==3)
            ratio_three = (dis[three_index,1])[0]
        except:
            ratio_three = 0

        ratio = ratio_final / (ratio_one+ratio_two+ratio_three)

    results['ratios'].append(ratio)

all_modes = np.array(results['mode'])
all_eis = np.array(results['incidence_es'])
all_ratios = np.array(results['ratios'])


for mode in ['ldfa','bomd','tdpt','pes']:
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
    if mode=='pes':
        mode_args = pes_args.copy()

    a = ax.plot(incidence_es,ratios,**mode_args,markersize=6,markeredgecolor='black')

# 
###########################

for i in range(2):
    for j in range(2):
        font='Arial'
        for tick in ax[i,j].get_xticklabels():
            tick.set_fontname(font)
        for tick in ax[i,j].get_yticklabels():
            tick.set_fontname(font)


        ax[i,j].tick_params(axis='both', which='major', labelsize=fontsize)
        ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i,j].xaxis.set_major_locator(MultipleLocator(0.2))



fig.set_figheight(3.)
fig.set_figwidth(3.25)

plt.gcf().subplots_adjust(left=0.3,bottom=0.3)
fig.savefig('fig3.pdf',transparent=True,bbox_inches='tight',dpi=300)
fig.savefig('fig3.eps',transparent=False)#,bbox_inches='tight')
