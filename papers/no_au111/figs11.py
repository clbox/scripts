from matplotlib.pyplot import xcorr
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
import seaborn as sns
import seaborn.distributions as sd
from seaborn.palettes import color_palette, blend_palette
from six import string_types

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

markers = ['o','x','s','v','.','x']
linestyles = ['-','--']

red = sns.color_palette("Reds")[-2]
blue = sns.color_palette("Blues")[-2]
colours = [blue,red,'darkgreen','goldenrod','violet','pink']


O_mass = 14.0067
N_mass = 15.999

marker_args = [{'linewidth' : 0.7, 's' : 10,'facecolor' : 'none','alpha' : 0.6},{'linewidth' : 0.7, 's' : 20, 'facecolor' : 'maroon', 'alpha' : 0.6}]

def _bivariate_kdeplot(x, y, filled, fill_lowest,
                       kernel, bw, gridsize, cut, clip,
                       axlabel, cbar, cbar_ax, cbar_kws, ax, **kwargs):
    """Plot a joint KDE estimate as a bivariate contour plot."""
    # Determine the clipping
    if clip is None:
        clip = [(-np.inf, np.inf), (-np.inf, np.inf)]
    elif np.ndim(clip) == 1:
        clip = [clip, clip]

    # Calculate the KDE
    if sd._has_statsmodels:
        xx, yy, z = sd._statsmodels_bivariate_kde(x, y, bw, gridsize, cut, clip)
    else:
        xx, yy, z = sd._scipy_bivariate_kde(x, y, bw, gridsize, cut, clip)

    # Plot the contours
    n_levels = kwargs.pop("n_levels", 10)
    cmap = kwargs.get("cmap", "BuGn" if filled else "BuGn_d")
    if isinstance(cmap, string_types):
        if cmap.endswith("_d"):
            pal = ["#333333"]
            pal.extend(color_palette(cmap.replace("_d", "_r"), 2))
            cmap = blend_palette(pal, as_cmap=True)
        else:
            cmap = plt.cm.get_cmap(cmap)

    kwargs["cmap"] = cmap
    contour_func = ax.contourf if filled else ax.contour
    cset = contour_func(xx, yy, z, n_levels, **kwargs)
    if filled and not fill_lowest:
        cset.collections[0].set_alpha(0)
    kwargs["n_levels"] = n_levels

    if cbar:
        cbar_kws = {} if cbar_kws is None else cbar_kws
        ax.figure.colorbar(cset, cbar_ax, ax, **cbar_kws)

    # Label the axes
    if hasattr(x, "name") and axlabel:
        ax.set_xlabel(x.name)
    if hasattr(y, "name") and axlabel:
        ax.set_ylabel(y.name)

    return ax, cset

mode2 = 'orient' #orient or final

results = {'final_v' : [],'atom_first' : [], 'pos1' : [], 'pos2' : []}
filenames = sys.argv[1:]
for i,filename in enumerate(filenames):
    trapped = False
    print(filename)
    ntrajs = 0
    v_f=-1
    with open(filename) as f:
        if 'trapped' in filename:
            trapped = True
        else:
            v_f = int(filename.replace('.dat',''))

        for line in f:

            if 'Trajectory' in line:
                ntrajs += 1
                
            elif 'Initial' in line:
                #NOT trapped
                numbers = line.replace(',',' ')
                i_v = float(numbers.split()[8])
                i_r = float(numbers.split()[9])
                i_t = float(numbers.split()[10])
            elif 'Final' in line:
                #NOT trapped
                numbers = line.replace(',',' ')
                f_v = float(numbers.split()[8])
                f_r = float(numbers.split()[9])
                f_t = float(numbers.split()[10])

            elif 'Impact 0:' in line:
                pos1 = [float(a) for a in line.split()[2:]]

            elif 'Impact 1:' in line:
                pos2 = [float(a) for a in line.split()[2:]]
            
            elif 'Lifetime' in line:
                #NOT trapped
                numbers = line.replace(',','')
                lifetime = float(numbers.split()[2])
                scat = float(numbers.split()[7])
                jf = int(numbers.split()[-1])

            elif 'first' in line:
                if 'N first' in  line:
                    atom_first = 'n'
                if 'O first' in  line:
                    atom_first = 'o'

            elif '-----------------' in line:
                if trapped:
                    results['final_v'].append(-1)
                else:
                    results['final_v'].append(v_f)
                
                results['atom_first'].append(atom_first)
                results['pos1'].append(pos1)
                results['pos2'].append(pos2)




weights = [1,0.2]
N_trajs_r = []
N_trajs_theta = []
cmaps = ['Blues','Reds']
kwargs = {'levels': np.arange(0, 0.15, 0.01)}
com_bins = np.linspace(1,3,100)
for mode in ['theta', 'r']:
    fig = plt.figure()

    gs = GridSpec(4,4)

    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_x = fig.add_subplot(gs[0,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])

    if mode2 == 'orient':
        labels = [r'N $\downarrow$',r'O $\downarrow$']
        for i,atom in enumerate(['n','o']):
            zorder=5-i
            indices = (np.where((np.array(results['final_v'])==-1) & (np.array(results['atom_first'])==atom)))[0]
            #indices = (np.where((np.array(results['atom_first'])==atom)))[0]
            print('Number of {} first is {}'.format(atom,len(indices)))
            #indices = (np.where((np.array(results['final_v'])==-1) & (np.array(results['atom_first'])==atom)))[0]
            pos1 = np.array((results['pos1']))[indices,:]
            pos2 = np.array((results['pos2']))[indices,:]

            y1 = pos1[:,2]
            y2 = pos2[:,2]
            centre_mass_z = (y1*O_mass + y2*N_mass) / (O_mass + N_mass)

            COM = (pos1*O_mass + pos2*N_mass) / (O_mass + N_mass)


            dx = pos1[:,0]-pos2[:,0]
            dy = pos1[:,1]-pos2[:,1]
            dz = pos1[:,2]-pos2[:,2]
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            angle = np.arcsin((y1-y2)/r) * 180/np.pi

            weights = np.ones_like(r)/len(r)
            sd._bivariate_kdeplot = _bivariate_kdeplot
            if mode == 'r':
                r_bins = np.linspace(1,2.,100)
                N_trajs_r.append(len(r))
                #ax_joint.scatter(r,centre_mass_z,zorder=zorder,label=labels[i],marker=markers[i], **marker_args[i],edgecolors=colours[i])
                _, cs = sns.kdeplot(r,centre_mass_z,cmap=cmaps[i],ax=ax_joint,shade=False,gridsize=100,cbar=False,shade_lowest=True,levels=5)#**kwargs)
                #sns.histplot(r,centre_mass_z,bins=50,ax=ax_joint)
                ax_marg_x.hist(r,bins=r_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
            else:
                theta_bins = np.linspace(-90,90,100)
                N_trajs_theta.append(len(r))
                _, cs  = sns.kdeplot(angle,centre_mass_z,cmap=cmaps[i],ax=ax_joint,shade=False,gridsize=100,cbar=False,shade_lowest=True,levels=5)#,**kwargs)
                #ax_joint.scatter(angle,centre_mass_z,zorder=zorder,label=labels[i],marker=markers[i], **marker_args[i],edgecolors=colours[i])
                ax_marg_x.hist(angle,bins=theta_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
            ax_marg_y.hist(centre_mass_z,bins=com_bins,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder,density=False)
            #plt.clabel(cs, cs.levels, inline=True,fontsize=4)
    else:
        print('nada')
            

    





    #Limits
    #ax_joint.legend(ncol=2,handletextpad=0.15,columnspacing=0.2,fancybox=True,framealpha=1,handlelength=1,bbox_to_anchor=(0.5, 1.05), loc='center')


    annotate_args = {'xy' : (0.01,0.05), 'xycoords' : 'axes fraction'}
    ax_joint.annotate(r'$v_i = 3$',ha="left", **annotate_args)
    annotate_args['xy'] = (0.01,0.9)
    if mode == 'r':
        ax_joint.set_xlim(1.,1.5)
        ax_marg_x.set_xlim(1.,1.5)
        ax_joint.annotate(r'(b)',ha="left", **annotate_args)
        ax_joint.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax_joint.xaxis.set_major_locator(MultipleLocator(0.1))
        ax_joint.set_xlabel(r"$r$ / $\mathrm{\AA{}}$")
        ax_marg_x.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax_marg_x.xaxis.set_major_locator(MultipleLocator(0.2))
        ax_marg_y.set_ylim(1.25,3)
        ax_joint.set_ylim(1.25,3)
        ax_joint.set_ylabel(r"COM height / $\mathrm{\AA{}}$")

        ax_marg_x.set_ylim(0,80)
        ax_marg_y.set_xlim(0,40)
        ax_marg_y.xaxis.set_major_locator(MultipleLocator(40))
        ax_marg_y.xaxis.set_minor_locator(MultipleLocator(10))

        ax_marg_x.yaxis.set_major_locator(MultipleLocator(80))
        ax_marg_x.yaxis.set_minor_locator(MultipleLocator(10))
        ax_joint.text(1.3, 1.5, "N$\downarrow$", size=12, color=blue)
        ax_joint.text(1.35, 2.8, "O$\downarrow$", size=12, color=red)

        ax_joint.text(1.3, 1.3, r"$n$ ="+str(N_trajs_r[0]), size=12, color=blue,ha='left')
        ax_joint.text(1.35, 2.6, r"$n$ ="+str(N_trajs_r[1]), size=12, color=red,ha='left')
        ax_marg_x.set_ylabel(r'Freq / $\mathrm{\AA{}}^{-1}$')
        ax_marg_y.set_xlabel(r'Freq / $\mathrm{\AA{}}^{-1}$')

    else:
        ax_marg_x.set_xlim(-90,90)
        ax_joint.set_xlim(-90,90)
        ax_joint.annotate(r'(a)',ha="left", **annotate_args)
        ax_joint.xaxis.set_minor_locator(MultipleLocator(5))
        ax_joint.xaxis.set_major_locator(MultipleLocator(30))
        ax_joint.set_xlabel(r"$\theta_\mathrm{{o}}$")
        ax_marg_x.xaxis.set_minor_locator(MultipleLocator(5))
        ax_marg_x.xaxis.set_major_locator(MultipleLocator(30))
        ax_marg_y.set_ylim(1.25,3)
        ax_joint.set_ylim(1.25,3)
        ax_joint.set_ylabel(r"COM height / $\mathrm{\AA{}}$")

        ax_marg_x.set_ylim(0,30)
        ax_marg_y.set_xlim(0,40)

        ax_marg_y.xaxis.set_major_locator(MultipleLocator(40))
        ax_marg_y.xaxis.set_minor_locator(MultipleLocator(10))
        ax_marg_x.yaxis.set_major_locator(MultipleLocator(30))
        ax_marg_x.yaxis.set_minor_locator(MultipleLocator(10))
        ax_joint.text(30, 1.5, "N$\downarrow$", size=12, color=blue)
        ax_joint.text(-30, 2.6, "O$\downarrow$", size=12, color=red)

        ax_joint.text(30, 1.3, r"$n$ ="+str(N_trajs_theta[0]), size=12, color=blue,ha='left')
        ax_joint.text(-30, 2.4, r"$n$ ="+str(N_trajs_theta[1]), size=12, color=red,ha='left')
        ax_marg_x.set_ylabel('Freq')
        ax_marg_y.set_xlabel(r'Freq / $\mathrm{\AA{}}^{-1}$')




    ax_joint.yaxis.set_major_locator(MultipleLocator(0.25))
    ax_joint.yaxis.set_minor_locator(MultipleLocator(0.05))

    ax_marg_y.yaxis.set_major_locator(MultipleLocator(0.25))
    ax_marg_y.yaxis.set_minor_locator(MultipleLocator(0.05))

    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    # Set labels on joint
    # 



    # 
    plt.subplots_adjust(hspace=0.6,wspace=0.6)
    fig.set_figwidth(3.25)
    fig.set_figheight(3.0)

    fig.savefig('figs11'+mode+'.pdf',transparent=True,bbox_inches='tight')
    fig.savefig('figs11'+mode+'.png',transparent=True,bbox_inches='tight',dpi=300)
    #fig.savefig('figs11'+mode+'.eps',transparent=False,bbox_inches='tight')
    fig.savefig('figs11'+mode+'scatter.tiff',dpi=600,transparent=True,bbox_inches='tight')