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
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

filenames = sys.argv[1:]


labels = ['Trapped',r'$\nu_f = 1$',r'$\nu_f = 2$',r'$\nu_f = 3$']

#labels = ['Trapped',r'Inelastic',r'Elastic']

O_heights = []
N_heights = []
Coms_z = []
bond_lengths = []
thetas = []

O_mass = 14.0067
N_mass = 15.999

markers = ['o','^','s','v']

colours = ['navy','maroon','darkgreen','goldenrod']
# Read the array from disk
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
for i,filename in enumerate(filenames):
    impact_geos = np.loadtxt(filename)

    # Note that this returned a 2D array!
    print(impact_geos.shape)

    impact_geos = impact_geos.reshape(impact_geos.shape[0]//2,2,3)

    print(impact_geos.shape)



    ############ With histograms ################ 
    x1 = impact_geos[:,0,0] #Ox
    y1 = impact_geos[:,0,2] #Oz

    x2 = impact_geos[:,1,0] #Nx
    y2 = impact_geos[:,1,2] #Nz

    O_heights.append(y1)
    N_heights.append(y2)
    

    centre_mass_z = (y1*O_mass + y2*N_mass) / (O_mass + N_mass)
    Coms_z.append(centre_mass_z)


    dx = impact_geos[:,0,0]-impact_geos[:,1,0]
    dy = impact_geos[:,0,1]-impact_geos[:,1,1]
    dz = impact_geos[:,0,2]-impact_geos[:,1,2]
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    bond_lengths.append(r)

    angle = np.arcsin((y1-y2)/r) * 180/np.pi
    thetas.append(angle)

    bins_x = np.linspace(-5, 7, 100)
    bins_y = np.linspace(0, 7, 100)

    # sns.distplot(y1,bins=bins_y,color='red')
    # sns.distplot(y1,bins=bins_y,color='blue')


    #ax.hist(y1,bins_y,color='red',alpha=0.5,label='Oxygen')
    #ax.hist(y2,bins_y,color='blue',alpha=0.5,label='Nitrogen')
    sns.distplot(y2, bins = bins_y, hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1}, 
                  label = labels[i], color = colours[i])

font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)


ax.tick_params(axis='both', which='major', labelsize=12)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.legend(fontsize=15)

ax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_major_locator(MultipleLocator(1.0))

ax.set_ylabel('Density',fontsize=12,fontname=font,color='black')
ax.set_xlabel(r"Nitrogen impact height / $\mathrm{\AA{}}$",fontsize=12,fontname=font)

fig.set_figheight(2.0)
fig.set_figwidth(3.25)
plt.gcf().subplots_adjust(left=0.3,bottom=0.3)

ax.get_legend().remove()
#ax.legend(loc=1)
ax.set_xlim(0,6.0)
fig.savefig('impact_histz_only.pdf',transparent=True,bbox_inches='tight')
#fig.savefig('_histz.png',dpi=300,transparent=True,bbox_inches='tight')


fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
for i in range(len(labels)):
    ax.scatter(O_heights[i],N_heights[i],marker=markers[i],label=labels[i],zorder=4-i,s=10,facecolors="None",edgecolors=colours[i])
ax.set_ylabel(r"Nitrogen impact height / $\mathrm{\AA{}}$",fontsize=12,fontname=font)
ax.set_xlabel(r"Oxygen impact height / $\mathrm{\AA{}}$",fontsize=12,fontname=font)

#ax.legend(loc=3,ncol=1,fontsize=12,fancybox=True,framealpha=0)

ax.set_xlim(0.5,4)
ax.set_ylim(0.5,4)




ax.tick_params(axis='both', which='major', labelsize=12)
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))

ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# ax.xaxis.set_major_locator(MultipleLocator(1.0))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
# ax.yaxis.set_major_locator(MultipleLocator(1.0))
ax.set_xticks([0.5,1.5,2.5,3.5])
ax.set_yticks([0.5,1.5,2.5,3.5])

fig.set_figheight(3.25)
fig.set_figwidth(3.25)
plt.gcf().subplots_adjust(left=0.3,bottom=0.3)

fig.savefig('heights_correlation.pdf',transparent=True,bbox_inches='tight')


fig.legend(ncol=4,fontsize=12,fancybox=True,framealpha=0)
fig.set_figwidth(5.25)
fig.set_figheight(1.7)
ax.remove()
fig.savefig('legend.pdf',transparent=True)#,bbox_inches='tight')


fig = plt.figure()

gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])

bins_x = np.linspace(-5, 7, 100)
bins_y = np.linspace(0, 7, 100)
for i in range(len(filenames)):
    if i == 3:
        zorder=6
    else:
        zorder=4-i
    ax_joint.scatter(O_heights[i],N_heights[i],zorder=zorder,s=10,label=labels[i],marker=markers[i],facecolors="None",edgecolors=colours[i])
    ax_marg_x.hist(O_heights[i],bins_x,color=colours[i],alpha=0.5,zorder=zorder)
    ax_marg_y.hist(N_heights[i],bins_y,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder)

#Limits
ax_marg_y.set_ylim(0.5,4)
ax_marg_x.set_xlim(0.5,4)
ax_joint.set_xlim(0.5,4)
ax_joint.set_ylim(0.5,4)

#ax_joint.legend(loc=0,ncol=1,fontsize=12,fancybox=True,framealpha=0)

# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel(r"Oxygen impact height / $\mathrm{\AA{}}$")
ax_joint.set_ylabel(r"Nitrogen impact height / $\mathrm{\AA{}}$")

# Set labels on marginals
ax_marg_y.set_xlabel('Frequency')
ax_marg_x.set_ylabel('Frequency')

fig.set_figwidth(3.25)
fig.set_figheight(3.0)

fig.savefig('height_corr_hist.pdf',transparent=True,bbox_inches='tight')




#COM bond length ########
fig = plt.figure()

gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])

bins_x = np.linspace(1, 1.5, 100)
bins_y = np.linspace(1.5, 3, 100)
for i in range(len(filenames)):
    if i == 3:
        zorder=6
    else:
        zorder=4-i
    ax_joint.scatter(bond_lengths[i],Coms_z[i],zorder=zorder,s=10,label=labels[i],marker=markers[i],facecolors="None",edgecolors=colours[i])
    ax_marg_x.hist(bond_lengths[i],bins_x,color=colours[i],alpha=0.5,zorder=zorder)
    ax_marg_y.hist(Coms_z[i],bins_y,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder)

#Limits
ax_marg_y.set_ylim(1.5,3)
ax_marg_x.set_xlim(1.1,1.45)
ax_joint.set_xlim(1.1,1.45)
ax_joint.set_ylim(1.5,3)

#ax_joint.legend(loc=0,ncol=1,fontsize=12,fancybox=True,framealpha=0)

# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel(r"Bond length / $\mathrm{\AA{}}$")
ax_joint.set_ylabel(r"COM height / $\mathrm{\AA{}}$")

# Set labels on marginals
ax_marg_y.set_xlabel('Frequency')
ax_marg_x.set_ylabel('Frequency')

fig.set_figwidth(3.25)
fig.set_figheight(3.0)

fig.savefig('com_r_corr_hist.pdf',transparent=True,bbox_inches='tight')




#COM theta ########
fig = plt.figure()

gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])

bins_x = np.linspace(-110, 110, 100)
bins_y = np.linspace(1.5, 3, 100)
for i in range(len(filenames)):
    if i == 3:
        zorder=6
    else:
        zorder=4-i
    ax_joint.scatter(thetas[i],Coms_z[i],zorder=zorder,s=10,label=labels[i],marker=markers[i],facecolors="None",edgecolors=colours[i])
    ax_marg_x.hist(thetas[i],bins_x,color=colours[i],alpha=0.5,zorder=zorder)
    ax_marg_y.hist(Coms_z[i],bins_y,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder)

#Limits
ax_marg_y.set_ylim(1.5,3)
ax_marg_x.set_xlim(-90,90)
ax_joint.set_xlim(-90,90)
ax_joint.set_ylim(1.5,3)



ax_joint.xaxis.set_minor_locator(MultipleLocator(5))
ax_joint.xaxis.set_major_locator(MultipleLocator(30))
ax_joint.yaxis.set_major_locator(MultipleLocator(0.25))
ax_joint.yaxis.set_minor_locator(MultipleLocator(0.05))

ax_marg_x.xaxis.set_minor_locator(MultipleLocator(5))
ax_marg_x.xaxis.set_major_locator(MultipleLocator(30))
ax_marg_y.yaxis.set_major_locator(MultipleLocator(0.25))
ax_marg_y.yaxis.set_minor_locator(MultipleLocator(0.05))


#ax_joint.legend(loc=0,ncol=1,fontsize=12,fancybox=True,framealpha=0)

# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel(r"$\theta$")
ax_joint.set_ylabel(r"COM height / $\mathrm{\AA{}}$")

# Set labels on marginals
ax_marg_y.set_xlabel('Frequency')
ax_marg_x.set_ylabel('Frequency')

fig.set_figwidth(3.25)
fig.set_figheight(3.0)

fig.savefig('com_theta_corr_hist.pdf',transparent=True,bbox_inches='tight')




# #3d
# from matplotlib import pyplot
# matplotlib.use( 'tkagg' )

# fig = pyplot.figure()
# ax = Axes3D(fig)

# for i in range(len(filenames)):
#     ax.scatter(thetas[i],bond_lengths[i],Coms_z[i],marker=markers[i],color=colours[i],label=labels[i])
# ax.legend()
# pyplot.show()