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
plt.style.use('clb_publication_thesis')
#COMMAND: python ~/Documents/scripts/papers/no_au111/fig2.py v02/translational/*/300K/640/states_e.txt {v03,v11,v16}/overview/*/states_e.txt && open fig2.pdf

tdpt_args = {'marker' : 'o', 'linestyle' : '--','color' : 'mediumorchid', 'label' : r'MDEF(ODF)', 'alpha' : 1.0}
tdpt_d4_args = {'marker' : '^', 'linestyle' : '-','color' : 'cyan', 'label' : r'MDEF(ODF) $\Lambda_{rr} \times 4$', 'alpha' : 1.0}
#tdpt_args = {'marker' : '^', 'linestyle' : '--','color' : 'grey', 'label' : r'ODF', 'alpha' : 1.0}
bomd_args = {'marker' : '^','linestyle' : '-','color' : 'red', 'label' : r'BOMD', 'alpha' : 1.0}
ldfa_args = {'zorder' : -2, 'marker' : 's','linestyle' : '-.','color' : 'blue', 'label' : r'MDEF(LDFA)', 'alpha' : 1.0}
ldfa_i4_args = {'zorder' : -1, 'marker' : 's','linestyle' : '-.','color' : 'pink', 'label' : r'MDEF(LDFA) $\mathbf{\Lambda} \times 4$', 'alpha' : 1.0}

mdef_args = {'zorder' : -4,'marker' : 's','linestyle' : '-','color' : '#F5C799', 'markerfacecolor' : 'white', 'label' : r'MDEF$^{[2]}$', 'alpha' : 1.0}
iesh_args = {'marker' : 'o','linestyle' : '-','color' : '#9ABD8F', 'markerfacecolor' : 'white', 'label' : r'IESH$^{[2]}$', 'alpha' : 1.0}
annotate_args = {'xy' : (0.93,0.85), 'xycoords' : 'axes fraction'}
annotate_args2 = {'xy' : (0.05,0.85), 'xycoords' : 'axes fraction'}
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

# fig, ax = plt.subplots(1, 2, sharex='all',sharey='all')#, constrained_layout=True)
fig, ax1 = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)
#fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)
#v02 exp
ax = [ax1]
ax[0].bar(x2_exp,v2_exp,color=exp_colour,edgecolor='black',label='Expt', zorder=-100)#label=r'$v_i=2$ exp')
ax[0].set_xlim(0,4)
ax[0].set_ylim(0,1.0)
ax[0].annotate(r'$v_i = 2$',ha="right", **annotate_args)


v2_bomd = np.loadtxt('v02/translational/bomd/300K/640/states_e.txt')
v2_ldfa = np.loadtxt('v02/translational/ldfa/300K/640/states_e.txt')
v2_odf = np.loadtxt('v02/translational/tdpt/300K/640/states_e.txt')

#a = ax[0].plot(v2_bomd[:,0],v2_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
c = ax[0].plot(v2_ldfa[:,0],v2_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)
b = ax[0].plot(v2_odf[:,0],v2_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)

ax[0].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[0].xaxis.set_major_locator(MultipleLocator(1))


ax[0].set_xlabel(r"$v_f$")
ax[0].set_ylabel('Population',)

handles,labels = ax[0].get_legend_handles_labels()
print(labels)
#plt.legend(ncol=2,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.1), loc='center')

fig.set_figheight(2.)
fig.set_figwidth(2.25)
fig.savefig('thesis.pdf',transparent=False,bbox_inches='tight')


######################### v3 #################################
########################################################################
########################################################################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
ax = [None,ax1]
#v03 exp
ax[1].bar(x3_exp,v3_exp,color=exp_colour,edgecolor='black',label=r'Expt')
ax[1].set_xlim(0,4)
ax[1].set_ylim(0,1.0)
ax[1].annotate(r'$v_i = 3$',ha="right", **annotate_args)


v3_bomd = np.loadtxt('v03/overview/bomd/states_e.txt')
v3_odf = np.loadtxt('v03/overview/tdpt/states_e.txt')
v3_ldfa = np.loadtxt('v03/overview/ldfa/states_e.txt')

a = ax[1].plot(v3_bomd[:,0],v3_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
c = ax[1].plot(v3_ldfa[:,0],v3_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)
b = ax[1].plot(v3_odf[:,0],v3_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)


ax[1].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[1].xaxis.set_major_locator(MultipleLocator(1))


ax[1].set_xlabel(r"$v_f$")
# ax[1].set_ylabel('Population',)

handles,labels = ax[1].get_legend_handles_labels()
print(labels)

plt.legend(ncol=4,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(-0.1, 1.1), loc='center')

fig.set_figheight(2.4)
fig.set_figwidth(5.25)
fig.savefig('thesis_2.pdf',transparent=False,bbox_inches='tight')

##########################################################################################
##########################################################################################
##########################################################################################


######################### v11 ############################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
#v02 exp
ax.bar(x11_exp,v11_exp,color=exp_colour,edgecolor='black',label=r'$v_i=11$ exp',zorder=-10)
ax.set_ylim(0,0.5)
ax.set_xlim(0,12)
ax.annotate(r'$v_i = 11$',ha="left", **annotate_args2)
ax.annotate(r'$E_i = 0.95\ \mathrm{eV}$',ha="left", xy=(0.05,0.75), xycoords='axes fraction')

v11_bomd = np.loadtxt('v11/overview/bomd/states_e.txt')
v11_odf = np.loadtxt('v11/overview/tdpt/states_e.txt')
v11_ldfa = np.loadtxt('v11/overview/ldfa/states_e.txt')
v11_iesh = np.loadtxt('v11/overview/ref/iesh/states_e.txt')
v11_mdef = np.loadtxt('v11/overview/ref/mdef/states_e.txt')
v11_odf_rr4 = np.loadtxt('v11/scaled/tdpt/d4/states_e.txt')

v11_ldfa_i4 = np.loadtxt('v11/scaled/ldfa/i4/states_e.txt')

a = ax.plot(v11_bomd[:,0],v11_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
b = ax.plot(v11_odf[:,0],v11_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)
c = ax.plot(v11_iesh[:,0],v11_iesh[:,1],markersize=6,markeredgecolor='black',**iesh_args)
d = ax.plot(v11_mdef[:,0],v11_mdef[:,1],markersize=6,markeredgecolor='black',**mdef_args)
e = ax.plot(v11_ldfa[:,0],v11_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)
f = ax.plot(v11_odf_rr4[:,0],v11_odf_rr4[:,1],markersize=6,markeredgecolor='black',**tdpt_d4_args)
g = ax.plot(v11_ldfa_i4[:,0],v11_ldfa_i4[:,1],markersize=6,markeredgecolor='black',**ldfa_i4_args)

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(1))


ax.set_xlabel(r"$v_f$")
ax.set_ylabel('Population',)

handles,labels = ax.get_legend_handles_labels()
print(labels)

#plt.legend(ncol=3,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.1), loc='center')

fig.set_figheight(3.)
fig.set_figwidth(4)
fig.savefig('thesis_3.pdf',transparent=False,bbox_inches='tight')


##########################################################################################
##########################################################################################

######################### v16 ############################

fig, ax = plt.subplots(1, 1)#, sharex='all',sharey='all')#, constrained_layout=True)
#v02 exp
ax.bar(x16_exp,v16_exp,color=exp_colour,edgecolor='black',label='Expt$^{[5]}$',zorder=-10)#,label=r'$v_i=16$ exp')
ax.set_ylim(0,0.25)
ax.set_xlim(0,18)
ax.annotate(r'$v_i = 16$',ha="left", **annotate_args2)
ax.annotate(r'$E_i = 0.52\ \mathrm{eV}$',ha="left", xy=(0.05,0.75), xycoords='axes fraction')

v16_bomd = np.loadtxt('v16/overview/bomd/states_e.txt')
v16_odf = np.loadtxt('v16/overview/tdpt/states_e.txt')
v16_iesh = np.loadtxt('v16/overview/ref/iesh/states_e.txt')
v16_mdef = np.loadtxt('v16/overview/ref/mdef/states_e.txt')
v16_ldfa = np.loadtxt('v16/overview/ldfa/states_e.txt')
v16_odf_rr4 = np.loadtxt('v16/scaled/tdpt/d4/states_e.txt')
v16_ldfa_i4 = np.loadtxt('v16/scaled/ldfa/i4/states_e.txt')

a = ax.plot(v16_bomd[:,0],v16_bomd[:,1],markersize=6,markeredgecolor='black',**bomd_args)
b = ax.plot(v16_odf[:,0],v16_odf[:,1],markersize=6,markeredgecolor='black',**tdpt_args)
c = ax.plot(v16_iesh[:,0],v16_iesh[:,1],markersize=6,markeredgecolor='black',**iesh_args)
d = ax.plot(v16_mdef[:,0],v16_mdef[:,1],markersize=6,markeredgecolor='black',**mdef_args)
e = ax.plot(v16_ldfa[:,0],v16_ldfa[:,1],markersize=6,markeredgecolor='black',**ldfa_args)
f = ax.plot(v16_odf_rr4[:,0],v16_odf_rr4[:,1],markersize=6,markeredgecolor='black',**tdpt_d4_args)
g = ax.plot(v16_ldfa_i4[:,0],v16_ldfa_i4[:,1],markersize=6,markeredgecolor='black',**ldfa_i4_args)

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(2))


ax.set_xlabel(r"$v_f$")
ax.set_ylabel('Population',)

handles,labels = ax.get_legend_handles_labels()
print(labels)
plt.legend(ncol=2,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.44, 1.2), loc='center')

fig.set_figheight(3.)
fig.set_figwidth(4)
fig.savefig('thesis_4.pdf',transparent=False,bbox_inches='tight')


##########################################################################################
##########################################################################################
