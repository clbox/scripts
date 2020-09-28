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
from mem_energy_loss import read_memory_kernel
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


x2_exp = np.arange(1,4,1)
v2_exp = [0.33,0.66,0.00038]

x3_exp = np.arange(1,4,1)
v3_exp = [0.2195652173913044,0.42391304347826086,0.35434782608695653]

x11_exp = np.arange(2,12)
v11_exp = [0.047,0.099,0.18,0.18,0.18,0.12,0.068,0.052,0.036,0.025]

x15_exp = np.arange(5,16)
v15_exp = [0.115,0.1339,0.194,0.192,0.125,0.082,0.04,0.05,0.019,0.015,0.036]

x16_exp = np.arange(0,17,1)
v16_exp = [0.0,0.0,0.04,0.08,0.13,0.15,0.19,0.11,0.12,0.07,0.04,0.02,0.03,0.02,0.01,0.02,0.02]

tdpt_args = {'marker' : 'o', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
#tdpt_args = {'marker' : '^', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'LDFA', 'alpha' : 1.0}

annotate_args = {'xy' : (0.96,0.80), 'xycoords' : 'axes fraction'}
exp_colour = 'gold'


#fig, ax = plt.subplots(2, 2, gridspec_kw={'width_ratios': [2,1]})#, sharex='all',sharey='all')#, constrained_layout=True)
fig = plt.figure()
gs0 = gridspec.GridSpec(3,1,height_ratios=[0.6,0.1,1])

gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0], hspace=0)
ax0 = fig.add_subplot(gs00[0])

gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[2], hspace=0.6)
ax1 = fig.add_subplot(gs01[0])
ax2 = fig.add_subplot(gs01[1])
# gs = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[1, 1, 1] )

# gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = gs[0],hspace=1)
# gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = gs[1],hspace=0.7)
# gs3 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = gs[2],hspace=0.7)




ax = np.array(([None,ax0],[ax1,ax2]))




#Plot spectrum
x0 = 0
s = 0.6
x = np.linspace(0,6,5000)
c1 = 'grey'
gauss = (np.exp(-0.5*((x-x0)*(x-x0))/(s*s))/(s*np.sqrt(np.pi)))*(1/np.sqrt(2))
ax[0,1].plot(x,gauss,'--',color=c1,linewidth=2)
# ax[0,1].fill_between(x, gauss,np.zeros_like(gauss),color=c1,alpha=0.4)

filename = 'spec/projected_memory_kernel.out'
bins,re,im,dimension,max_e = read_memory_kernel(filename,treat_complex=False)



output_dir = os.path.dirname(filename)
if 'project' in filename:
    tensor_file = (glob.glob(output_dir+'/*projected_tensor.out'))[0]

tensor = np.loadtxt(tensor_file)
element_val = tensor[0,0]

ax[0,1].annotate('', xy=(0,element_val), xycoords='data', xytext=(-0.3, element_val), 
        arrowprops=dict(arrowstyle="-|>, head_width=0.3, head_length=0.7",color=c1),
        #arrowprops=dict(width=0.5),
        color='red')

c2 = 'mediumorchid'
ax[0,1].axhline(y=(element_val*4), xmin=0, xmax=100,color=c2,linestyle=':',linewidth=1.5)

# ax[0,1].annotate('', xy=(2,4*element_val), xycoords='data', xytext=(2,element_val), 
#         arrowprops=dict(arrowstyle="-|>, head_width=0.3, head_length=0.7",color=c2),
#         #arrowprops=dict(width=0.5),
#         color='red')
ax[0,1].plot(bins,re[0,:],linestyle='-',linewidth=1,color='black')
ax[0,1].text( x=1.8,y=0.85, s=r'$\times 4$', color=c2)
ax[0,1].set_ylim(bottom=0,top=1.2)
ax[0,1].set_xlim(0,np.max(bins))
ax[0,1].xaxis.set_major_locator(MultipleLocator(0.3))
ax[0,1].yaxis.set_major_locator(MultipleLocator(0.3))
ax[0,1].xaxis.set_minor_locator(MultipleLocator(0.1))
ax[0,1].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[0,1].set_ylabel(r'$\Lambda_{rr}(\epsilon)$ / ps$^{-1}$',color='black')
ax[0,1].set_xlabel(r'Excitation energy / eV',color='black')


#fig.delaxes(ax[0,0])
#fig.delaxes(ax[0,1])
#v11
ax[1,0].bar(x11_exp,v11_exp,color=exp_colour,edgecolor='black',label=r'$v_i=11$ exp')
ax[1,0].set_ylim(0,0.3)
plotted_exp = True
ax[1,0].set_xlim(0,12)
ax[1,0].annotate(r'$v_i = 11$',ha="right", **annotate_args)

#v16
ax[1,1].bar(x16_exp,v16_exp,color=exp_colour,edgecolor='black',label='EXPT')#,label=r'$v_i=16$ exp')
ax[1,1].set_ylim(0,0.3)
ax[1,1].set_xlim(0,18)
ax[1,1].annotate(r'$v_i = 16$',ha="right", **annotate_args)



annotate_args['xy'] = (0.01,0.85)
# ax[0,0].annotate(r'(a)',ha="left", **annotate_args)
ax[0,1].annotate(r'(a)',ha="left", **annotate_args)
annotate_args['xy'] = (0.01,0.80)
ax[1,0].annotate(r'(b)',ha="left", **annotate_args)
ax[1,1].annotate(r'(c)',ha="left", **annotate_args)
indices=[]
for i,filename in enumerate(filenames):

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
        mode_args['label'] = mode_args['label'] + r'[$ \mathbf{\Lambda} \times 2$]'
        mode_args['linestyle'] = '--'

    if 'i3' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'[$ \mathbf{\Lambda} \times 3$]'
        mode_args['linestyle'] = ':'
    
    if 'i4' in os.path.abspath(filename):
        mode_args['label'] = mode_args['label'] + r'[$ \mathbf{\Lambda} \times 4$]'
        

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
        mode_args['label'] = mode_args['label'] + r'[$\Lambda_{rr} \times 4$]'
        mode_args['linestyle'] = '-'
        mode_args['marker'] = '^'
        mode_args['color'] = 'mediumorchid'

    if 'pes' in os.path.abspath(filename):
        mode_args['color'] = 'green'
        mode_args['marker'] = 'v'
        mode_args['label'] = mode_args['label'] + r'+ PES(2)'
    
    if 'v02' in os.path.abspath(filename):
        indices = [0,0]
    if 'v03' in os.path.abspath(filename):
        indices = [0,1]
    if 'v11' in os.path.abspath(filename):
        indices = [1,0]
    if 'v16' in os.path.abspath(filename):
        indices = [1,1]

    a = ax[indices[0],indices[1]].plot(dis[:,0],dis[:,1],**mode_args,markersize=4,markeredgecolor='black')


font='Arial'
for i in [1]:
    
    #ax[i,1].yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
    for j in range(2):
        ax[i,j].tick_params(axis='both', which='major')
        ax[i,j].xaxis.set_major_locator(MaxNLocator(integer=True))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(1))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.025))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.1))
        ax[i,j].set_xlabel(r"$v_f$",fontname=font)

        for tick in ax[i,j].get_xticklabels():
            tick.set_fontname(font)
        for tick in ax[i,j].get_yticklabels():
            tick.set_fontname(font)

ax[1,0].set_ylabel('Population',fontname=font,color='black')
ax[1,1].set_ylabel('Population',fontname=font,color='black')
ax[1,1].xaxis.set_major_locator(MultipleLocator(3))
fig.set_figheight(5)
fig.set_figwidth(3.25)
plt.legend(ncol=4,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 2.8), loc='center')
#plt.tight_layout()
#plt.subplots_adjust(hspace=1.2)
#plt.gcf().subplots_adjust(right=0.01)
fig.savefig('fig5.pdf',transparent=True,bbox_inches='tight')
fig.savefig('fig5.tiff',transparent=True,bbox_inches='tight',dpi=600)
fig.savefig('fig5.eps',transparent=True,bbox_inches='tight')
