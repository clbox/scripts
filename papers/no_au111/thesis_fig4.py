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

plt.style.use('clb_publication_thesis')


# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
annotate=True


#First argument, a = all results, e = exp only normalisation
# new limits for some plots


filenames = sys.argv[1:]
tdpt_args = {'marker' : 'o', 'linestyle' : '--','color' : 'mediumorchid', 'label' : r'MDEF(ODF)', 'alpha' : 1.0}
#tdpt_args = {'marker' : '^', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'MDEF(LDFA)', 'alpha' : 1.0}

annotate_args = {'xy' : (0.05,0.65), 'xycoords' : 'axes fraction'}
exp_colour = 'gold'

if os.path.exists("figs8.txt"):
    os.remove("figs8.txt")
fig, axes = plt.subplots(3, 1, sharex='all',sharey='all')#, constrained_layout=True)

ax = np.array(([None,axes[0]],[axes[1],axes[2]]))

#fig.delaxes(ax[0,0])
#ax[0,0].axis('off')
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
#ax[0,0].annotate(r'(a)',ha="left", **annotate_args)
ax[0,1].annotate(r'(a)',ha="left", **annotate_args)
ax[1,0].annotate(r'(b)',ha="left", **annotate_args)
ax[1,1].annotate(r'(c)',ha="left", **annotate_args)
indices=[]
for i,filename in enumerate(filenames):
    mode =''
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
        mode_args['label'] = mode_args['label'] + r'[$ \mathbf{\Lambda} \times 2$])'
        mode_args['linestyle'] = ':'
        mode_args['color'] = 'indigo'
        mode_args['marker'] = 'x'
        if 'tdpt' in os.path.abspath(filename):
            mode='ODF*2'

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
        mode_args['color'] = 'orange'
        mode_args['marker'] = 'D'
        mode_args['linestyle'] = '-.'
        mode_args['label'] = mode_args['label'] + r')[RS]'
        if 'tdpt' in os.path.abspath(filename):
            mode = 'ODF[RS]'
    
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
        t = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markeredgecolor='black')
    elif 'ldfa' in os.path.abspath(filename):
        l = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markeredgecolor='black')
    elif 'bomd' in os.path.abspath(filename):
        b = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markeredgecolor='black')

    if mode in ['BOMD','LDFA','ODF','ODF[RS]','ODF*2']:
        with open('figs8.txt','a+') as f:
            f.write('initial vib '+'3 ' + mode2 +'\n')
            f.write(mode + '\n')
            for s in range(np.shape(dis)[0]):
                f.write(str(int(dis[s,0]))+'    ')
                f.write(str(dis[s,1]))
                f.write('\n')


#ax[1,1].yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

for i in range(2):
    
    
    for j in range(2):
        if i == 0 and j==0:
            continue
        ax[i,j].tick_params(axis='both', which='major')
        ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(1))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.1))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.3))

        
ax[1,1].set_xlabel('$v_\mathrm{f}$',color='black')    
ax[0,1].set_ylabel('Population',color='black')
ax[1,0].set_ylabel('Population',color='black')
ax[1,1].set_ylabel('Population',color='black')

#ax[0,1].set_ylabel('Population',fontname=font,color='black')
# fig.text(0.5, 0.00, r"$v_f$", ha='center')
# ax[1,1].xaxis.set_major_locator(MultipleLocator(1))
# ax[0,1].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
# ax[1,0].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())



handles,labels = ax[0,1].get_legend_handles_labels()
print(labels)
# order = np.array([0,3,
#         1,4,
#         2])

# handles = [handles[i] for i in order]
# labels = [labels[i] for i in order]

# handles.append((p1,p2,p3))
# labels.append('Exp')

# labels[2] = labels[2]+')'
# labels[4] = labels[4]+')'
# plt.legend(handles,labels,numpoints=1,
#                 handler_map={tuple: HandlerTuple(ndivide=None)},
#                 ncol=2,handletextpad=0.15,columnspacing=0.6,
#                 fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.4, 3.5), loc='center')

# ax[0,1].yaxis.set_label_position("right")
# ax[0,1].yaxis.tick_right()

ax[0,1].legend(handles,labels,numpoints=1,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                ncol=2,handletextpad=0.15,columnspacing=0.6,
                fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.4, 1.2), loc='center')

# plt.legend(ncol=4,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(-0.2, 2.25), loc='center')
#plt.tight_layout()
plt.subplots_adjust(hspace=0.12,wspace=0.1)
#plt.gcf().subplots_adjust(right=0.01)
fig.set_figheight(5)
fig.set_figwidth(1.5748)
fig.savefig('thesis_fig4.pdf',transparent=True,bbox_inches='tight')
# fig.savefig('fig4_all.tiff',transparent=True,bbox_inches='tight',dpi=600)
# fig.savefig('fig4_all.eps',transparent=True,bbox_inches='tight')
