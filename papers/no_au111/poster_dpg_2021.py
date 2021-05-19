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
#COMMAND: python ~/Documents/scripts/papers/no_au111/fig2.py v02/translational/*/300K/640/states_1_e.txt {v03,v11,v16}/overview/*/states_1_e.txt && open fig2.pdf
SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#plt.style.use('dark_background')
# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
annotate=True
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"

tdpt_args = {'marker' : 'o', 'linestyle' : '--','color' : 'mediumorchid', 'label' : r'MDEF(ODF)', 'alpha' : 1.0}
#tdpt_args = {'marker' : '^', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'MDEF(LDFA)', 'alpha' : 1.0}
mdef_args = {'marker' : 's','linestyle' : '-','color' : '#F5C799', 'markerfacecolor' : 'white', 'label' : r'MDEF$^{2}$', 'alpha' : 1.0}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : '#9ABD8F', 'markerfacecolor' : 'white', 'label' : r'IESH$^{2}$', 'alpha' : 1.0}
annotate_args = {'xy' : (0.98,0.8), 'xycoords' : 'axes fraction'}
exp_colour = 'gold'

#################### EXP ####################################
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

######################### v02 ############################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
#v02 exp
ax.bar(x2_exp,v2_exp,color=exp_colour,edgecolor='black',label='Expt')#label=r'$v_i=2$ exp')
ax.set_xlim(0,4)
ax.set_ylim(0,1.0)
ax.annotate(r'$v_i = 2$',ha="right", **annotate_args)


v2_bomd = np.loadtxt('v02/translational/bomd/300K/640/states_1_e.txt')
v2_ldfa = np.loadtxt('v02/translational/ldfa/300K/640/states_1_e.txt')
v2_odf = np.loadtxt('v02/translational/tdpt/300K/640/states_1_e.txt')

a = ax.plot(v2_bomd[:,0],v2_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
c = ax.plot(v2_ldfa[:,0],v2_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)
b = ax.plot(v2_odf[:,0],v2_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(1))


ax.set_xlabel(r"$v_f$")
ax.set_ylabel('Population',)

handles,labels = ax.get_legend_handles_labels()
print(labels)
#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.1), loc='center')

fig.set_figheight(2.)
fig.set_figwidth(2.25)
fig.savefig('poster_dpg_2021.pdf',transparent=True,bbox_inches='tight')


######################### v3 #################################
########################################################################
########################################################################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
#v03 exp
ax.bar(x3_exp,v3_exp,color=exp_colour,edgecolor='black',label=r'$v_i=3$ exp')
ax.set_xlim(0,4)
ax.set_ylim(0,1.0)
ax.annotate(r'$v_i = 3$',ha="right", **annotate_args)


v3_bomd = np.loadtxt('v03/overview/bomd/states_1_e.txt')
v3_odf = np.loadtxt('v03/overview/tdpt/states_1_e.txt')
v3_ldfa = np.loadtxt('v03/overview/ldfa/states_1_e.txt')

a = ax.plot(v3_bomd[:,0],v3_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
b = ax.plot(v3_odf[:,0],v3_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)
c = ax.plot(v3_ldfa[:,0],v3_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(1))


ax.set_xlabel(r"$v_f$")
ax.set_ylabel('Population',)

handles,labels = ax.get_legend_handles_labels()
print(labels)
#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.1), loc='center')

fig.set_figheight(2.)
fig.set_figwidth(2.25)
fig.savefig('poster_dpg_2021_2.pdf',transparent=True,bbox_inches='tight')

##########################################################################################
##########################################################################################
##########################################################################################


######################### v11 ############################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
#v02 exp
ax.bar(x11_exp,v11_exp,color=exp_colour,edgecolor='black',label=r'$v_i=11$ exp')
ax.set_ylim(0,0.6)
ax.set_xlim(0,12)
ax.annotate(r'$v_i = 11$',ha="right", **annotate_args)


v11_bomd = np.loadtxt('v11/overview/bomd/states_1_e.txt')
v11_odf = np.loadtxt('v11/overview/tdpt/states_1_e.txt')
v11_ldfa = np.loadtxt('v11/overview/ldfa/states_1_e.txt')
v11_iesh = np.loadtxt('v11/overview/ref/iesh/states_1_e.txt')
v11_mdef = np.loadtxt('v11/overview/ref/mdef/states_1_e.txt')

a = ax.plot(v11_bomd[:,0],v11_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
b = ax.plot(v11_odf[:,0],v11_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)
c = ax.plot(v11_iesh[:,0],v11_iesh[:,1],markersize=6,markeredgecolor='black',**iesh_args)
d = ax.plot(v11_mdef[:,0],v11_mdef[:,1],markersize=6,markeredgecolor='black',**mdef_args)
#e = ax.plot(v11_ldfa[:,0],v11_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(1))


ax.set_xlabel(r"$v_f$")
ax.set_ylabel('Population',)

handles,labels = ax.get_legend_handles_labels()
print(labels)
#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.1), loc='center')

fig.set_figheight(2.5)
fig.set_figwidth(4.25)
fig.savefig('poster_dpg_2021_3.pdf',transparent=True,bbox_inches='tight')


##########################################################################################
##########################################################################################

######################### v16 ############################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
#v02 exp
ax.bar(x16_exp,v16_exp,color=exp_colour,edgecolor='black',label='Expt')#,label=r'$v_i=16$ exp')
ax.set_ylim(0,0.4)
ax.set_xlim(0,18)
ax.annotate(r'$v_i = 16$',ha="right", **annotate_args)


v16_bomd = np.loadtxt('v16/overview/bomd/states_1_e.txt')
v16_odf = np.loadtxt('v16/overview/tdpt/states_1_e.txt')
v16_iesh = np.loadtxt('v16/overview/ref/iesh/states_1_e.txt')
v16_mdef = np.loadtxt('v16/overview/ref/mdef/states_1_e.txt')
v16_ldfa = np.loadtxt('v16/overview/ldfa/states_1_e.txt')

a = ax.plot(v16_bomd[:,0],v16_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
b = ax.plot(v16_odf[:,0],v16_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)
c = ax.plot(v16_iesh[:,0],v16_iesh[:,1],markersize=6,markeredgecolor='black',**iesh_args)
d = ax.plot(v16_mdef[:,0],v16_mdef[:,1],markersize=6,markeredgecolor='black',**mdef_args)
#e = ax.plot(v16_ldfa[:,0],v16_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(2))


ax.set_xlabel(r"$v_f$")
ax.set_ylabel('Population',)

handles,labels = ax.get_legend_handles_labels()
print(labels)
plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.4), loc='center')

fig.set_figheight(2.5)
fig.set_figwidth(4.25)
fig.savefig('poster_dpg_2021_4.pdf',transparent=True,bbox_inches='tight')


##########################################################################################
##########################################################################################
