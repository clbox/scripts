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
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.image as mpimg
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
#COMMAND: python ~/Documents/scripts/papers/no_au111/fig2.py v02/translational/*/300K/640/states_1_e.txt {v03,v11,v16}/overview/*/states_1_e.txt && open fig2.pdf
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

#First argument, a = all results, e = exp only normalisation
# new limits for some plots


filenames = sys.argv[1:]
tdpt_args = {'marker' : 'o', 'linestyle' : '--','color' : 'mediumorchid', 'label' : r'ODF', 'alpha' : 1.0}
#tdpt_args = {'marker' : '^', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}

annotate_args = {'xy' : (0.05,0.65), 'xycoords' : 'axes fraction'}
exp_colour = 'gold'

if os.path.exists("fig4.txt"):
    os.remove("fig4.txt")
fig, ax = plt.subplots(2, 2)#, sharex='all',sharey='all')#, constrained_layout=True)


#fig.delaxes(ax[0,0])
ax[0,0].axis('off')
# img = mpimg.imread('drawing.png')
# ax[0,0].imshow(img)


#v03 - ISO
exp = np.loadtxt('v03_iso_950.txt')
err = np.array([0.189,0.206,0.473,0.306])-exp[:,1]
p1 = ax[0,1].bar(exp[:,0],exp[:,1],color=exp_colour,edgecolor='black'
        ,yerr=err,capsize=3, error_kw={'elinewidth' : 1})#,label=r'Expt')
ax[0,1].set_xlim(-0.5,3.5)
ax[0,1].set_ylim(0,0.9)
ax[0,1].annotate(r'Isotropic',ha="left",**annotate_args)

#v03 - N frst
exp = np.loadtxt('v03_nfirst_950.txt')
err = np.array([0.232,0.262,0.646,0.180])-exp[:,1]
p2 = ax[1,0].bar(exp[:,0],exp[:,1],color='cornflowerblue',edgecolor='black',
        yerr=err,capsize=3, error_kw={'elinewidth' : 1})#,label=r'Expt')
ax[1,0].set_xlim(-0.5,3.5)
ax[1,0].set_ylim(0,0.9)
ax[1,0].annotate(r'N$\downarrow$', ha="left",**annotate_args)

#v03 - O first
exp = np.loadtxt('v03_ofirst_950.txt')
err = np.array([0.143,0.165,0.401,0.456])-exp[:,1]
p3 = ax[1,1].bar(exp[:,0],exp[:,1],color="lightcoral",edgecolor='black',
        yerr=err,capsize=3, error_kw={'elinewidth' : 1})#,label=r'Expt'))#,label=r'Expt')
ax[1,1].set_xlim(-0.5,3.5)
ax[1,1].set_ylim(0,0.9)
ax[1,1].annotate(r'O$\downarrow$', ha="left",**annotate_args)


annotate_args['xy'] = (0.05,0.85)
ax[0,0].annotate(r'(a)',ha="left", **annotate_args)
ax[0,1].annotate(r'(b)',ha="left", **annotate_args)
ax[1,0].annotate(r'(c)',ha="left", **annotate_args)
ax[1,1].annotate(r'(d)',ha="left", **annotate_args)
indices=[]
for i,filename in enumerate(filenames):
    mode = ''
    dis = np.loadtxt(filename)
    mode_args = None
    if 'tdpt' in os.path.abspath(filename):
        mode_args = tdpt_args.copy()
        mode = 'ODF'

    if 'bomd' in os.path.abspath(filename):
        mode_args = bomd_args.copy()
        mode = 'BOMD'

    if 'ldfa' in os.path.abspath(filename):
        mode_args = ldfa_args.copy()
        mode = 'LDFA'

    # if '_1' in os.path.abspath(filename):
    #     mode_args['label'] = mode_args['label'] + ' SB'
    #     mode_args['linestyle'] = '--'
    #     mode_args['marker'] = 'x'

    if '_2' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + ' DB'
        mode_args['linestyle'] = ':'
        mode_args['marker'] = 'v'

    if 'i2' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'[$ \mathbf{\Lambda} \times 2$]'
        mode_args['linestyle'] = '-'
        mode_args['color'] = 'indigo'
        mode_args['marker'] = '^'

    if 'i3' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 3$'
        mode_args['linestyle'] = ':'
    
    if 'i4' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 4$'
        

        if 'tdpt' in os.path.abspath(filename):
            mode_args['color'] = 'mediumorchid'
            mode_args['marker'] = 'o'
            mode_args['linestyle'] = '-'
        if 'ldfa' in os.path.abspath(filename):
            mode_args['color'] = 'dodgerblue'
            mode_args['linestyle'] = '-.'

    if '_multi' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + ' MB'
        mode_args['linestyle'] = '-.'

    if 'd4' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r' (r) $\times 4$'
        mode_args['linestyle'] = ':'
        mode_args['marker'] = '^'
    
    if 'pes' in os.path.abspath(filename):
        mode_args['color'] = 'green'
        mode_args['marker'] = 'v'
        mode_args['linestyle'] = '-.'
        mode_args['label'] = mode_args['label'] + r' PES$_{\mathrm{rs}}$'
    
    if '_n' in os.path.abspath(filename):
        indices = [1,0]
        mode2 = 'N first'
    elif '_o' in os.path.abspath(filename):
        indices = [1,1]
        mode2 = 'O first'
    else:
        indices = [0,1]
        mode2 = 'Isotropic'


    if 'tdpt' in os.path.abspath(filename):
        t = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markersize=4,markeredgecolor='black')
    elif 'ldfa' in os.path.abspath(filename):
        l = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markersize=4,markeredgecolor='black')
    elif 'bomd' in os.path.abspath(filename):
        b = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markersize=4,markeredgecolor='black')

    if mode in ['BOMD','LDFA','ODF']:
        with open('fig4.txt','a+') as f:
            f.write('initial vib '+'3 ' + mode2 +'\n')
            f.write(mode + '\n')
            for s in range(np.shape(dis)[0]):
                f.write(str(int(dis[s,0]))+'    ')
                f.write(str(dis[s,1]))
                f.write('\n')


ax[1,1].yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
for i in range(2):
    ax[1,0].set_ylabel('Population',color='black')
    
    for j in range(2):
        if i == 0 and j==0:
            continue
        ax[i,j].tick_params(axis='both', which='major')
        ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(1))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.1))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.3))
        

        ax[0,j].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())


ax[0,1].set_ylabel('Population',color='black')
fig.text(0.5, 0.00, r"$v_f$", ha='center')
ax[1,1].xaxis.set_major_locator(MultipleLocator(1))




handles,labels = ax[0,1].get_legend_handles_labels()

handles.append((p1,p2,p3))
labels.append('Expt')

ax[0,1].legend(handles,labels,numpoints=1,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                ncol=4,handletextpad=0.15,columnspacing=0.6,
                fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(-0.0, 1.1), loc='center')

ax[0,1].yaxis.set_label_position("right")
ax[0,1].yaxis.tick_right()

fig.set_figheight(3.)
fig.set_figwidth(3.25)
# plt.legend(ncol=4,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(-0.2, 2.25), loc='center')
#plt.tight_layout()
plt.subplots_adjust(hspace=0.1,wspace=0.1)
#plt.gcf().subplots_adjust(right=0.01)
fig.savefig('fig4.pdf',transparent=True,bbox_inches='tight')
fig.savefig('fig4.tiff',transparent=True,bbox_inches='tight',dpi=600)
fig.savefig('fig4.eps',transparent=True,bbox_inches='tight')
