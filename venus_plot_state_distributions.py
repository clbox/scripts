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

x2_exp = np.arange(1,4,1)
v2_exp = [0.33,0.66,0.00038]

x3_exp = np.arange(0,7,1)
v3_exp = [0.0,0.2195652173913044,0.42391304347826086,0.35434782608695653,0.0,0,0]

x11_exp = np.arange(2,12)
v11_exp = [0.047,0.099,0.18,0.18,0.18,0.12,0.068,0.052,0.036,0.025]

x15_exp = np.arange(5,16)
v15_exp = [0.115,0.1339,0.194,0.192,0.125,0.082,0.04,0.05,0.019,0.015,0.036]

x16_exp = np.arange(0,17,1)
v16_exp = [0.0,0.0,0.04,0.08,0.13,0.15,0.19,0.11,0.12,0.07,0.04,0.02,0.03,0.02,0.01,0.02,0.02]

#tdpt_args = {'marker' : 'o', 'linestyle' : '-','color' : 'purple', 'label' : r'ODF', 'alpha' : 1.0}
tdpt_args = {'marker' : '^', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}

annotate_args = {'fontsize' : 12, 'xy' : (0.7,0.8), 'xycoords' : 'figure fraction'}
exp_colour = 'gold'


plotted_exp = False
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)
for i,filename in enumerate(filenames):

    if 'v02' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x2_exp,v2_exp,color=exp_colour,edgecolor='black',label='Exp')#label=r'$\nu_i=2$ exp')
        ax.set_xlim(0,4)
        ax.set_ylim(0,1.0)
        plotted_exp = True
        if annotate:
            ax.annotate(r'$\nu_i = 2$', **annotate_args)
    
    if 'v03' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x3_exp,v3_exp,color=exp_colour,edgecolor='black',label=r'$\nu_i=3$ exp')
        ax.set_xlim(0,6)
        ax.set_ylim(0,1.0)
        plotted_exp = True
        if annotate:
            ax.annotate(r'$\nu_i = 3$', **annotate_args)

    if 'v11' in os.path.abspath(filename) and not plotted_exp:
        #ax.bar(x11_exp,v11_exp,color='black',label=r'$\nu_i=11$ exp')
        ax.bar(x11_exp,v11_exp,color=exp_colour,edgecolor='black',label=r'$\nu_i=11$ exp')
        ax.set_xlim(0,16)
        ax.set_ylim(0,0.3)
        plotted_exp = True
        if annotate:
            ax.annotate(r'$\nu_i = 11$', **annotate_args)

    if 'v15' in os.path.abspath(filename) and not plotted_exp:
        ax.bar(x15_exp,v15_exp,color=exp_colour,edgecolor='black',label=r'$\nu_i=15$ exp')
        ax.set_xlim(0,17)
        ax.set_ylim(0,0.8)
        plotted_exp = True
        if annotate:
            ax.annotate(r'$\nu_i = 15$', **annotate_args)

    if 'v16' in os.path.abspath(filename) and not plotted_exp:
        #ax.bar(x16_exp,v16_exp,color='black',label=r'$\nu_i=16$ exp')
        ax.bar(x16_exp,v16_exp,color=exp_colour,edgecolor='black',label='Exp')#,label=r'$\nu_i=16$ exp')
        ax.set_ylim(0,0.3)
        ax.set_xlim(0,20)
        plotted_exp = True
        if annotate:
            ax.annotate(r'$\nu_i = 16$', **annotate_args)
        #ax.yaxis.set_major_formatter(plt.NullFormatter())

    dis = np.loadtxt(filename)
    mode_args = None
    if 'tdpt' in os.path.abspath(filename):
        mode_args = tdpt_args.copy()

    if 'bomd' in os.path.abspath(filename):
        mode_args = bomd_args.copy()

    if 'ldfa' in os.path.abspath(filename):
        mode_args = ldfa_args.copy()

    # if '_1' in os.path.abspath(filename):
    #     mode_args['label'] = mode_args['label'] + ' SB'
    #     mode_args['linestyle'] = '--'
    #     mode_args['marker'] = 'x'

    if '_2' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + ' DB'
        mode_args['linestyle'] = ':'
        mode_args['marker'] = 'v'

    if 'i2' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 2$'
        mode_args['linestyle'] = '--'

    if 'i3' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 3$'
        mode_args['linestyle'] = ':'
    
    if 'i4' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'$\times 4$'
        mode_args['linestyle'] = '-'

        if 'tdpt' in os.path.abspath(filename):
            mode_args['color'] = 'mediumorchid'
            mode_args['marker'] = 'o'
        if 'ldfa' in os.path.abspath(filename):
            mode_args['color'] = 'dodgerblue'

    if '_multi' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + ' MB'
        mode_args['linestyle'] = '-.'

    if 'd4' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r' (r) $\times 4$'
        mode_args['linestyle'] = ':'
        mode_args['marker'] = '^'
    
    if 'pes' in os.path.abspath(filename):
        mode_args['color'] = 'grey'


    

    a = ax.plot(dis[:,0],dis[:,1],**mode_args,markersize=6,markeredgecolor='black')

font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)


ax.tick_params(axis='both', which='major', labelsize=12)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=15)

ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.025))

fig.set_figheight(2.7)
fig.set_figwidth(3.25)
#fig.set_constrained_layout_pads(w_pad=0, h_pad=0)
ax.set_xlabel(r"Final vibrational state ($\nu_f$)",fontsize=12,fontname=font)
#ax.set_ylabel('Population',fontsize=12,fontname=font)#,color='white')
ax.set_ylabel('Population',fontsize=12,fontname=font,color='black')
plt.gcf().subplots_adjust(left=0.2,bottom=0.2)
#fig.text(0.5, 0.00, r"Final vibrational state ($\nu_f$)", ha='center',fontsize=15)
#fig.text(0.01, 0.5, 'Population', va='center', rotation='vertical',fontsize=15)
fig.savefig('state2state.pdf',transparent=True)#,bbox_inches='tight')


fig.legend(ncol=2,fontsize=12,fancybox=True,framealpha=0)
fig.set_figwidth(3.25)
fig.set_figheight(1.7)
ax.remove()
fig.savefig('legend.pdf',transparent=True)#,bbox_inches='tight')