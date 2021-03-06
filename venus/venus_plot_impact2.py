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

markers = ['o','^','s','v','.','x']
colours = ['navy','maroon','darkgreen','goldenrod','violet','pink']

O_mass = 14.0067
N_mass = 15.999


mode = sys.argv[1] #theta or r or xy

mode2 = sys.argv[2] #orient or final

results = {'final_v' : [],'atom_first' : [], 'pos1' : [], 'pos2' : []}
filenames = sys.argv[3:]
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



fig = plt.figure()

gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])




com_bins = np.linspace(1,3,100)

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

        if mode == 'r':
            r_bins = np.linspace(1,2.,100)
            ax_joint.scatter(r,centre_mass_z,zorder=zorder,s=1,label=labels[i],marker=markers[i],facecolors="None",edgecolors=colours[i])
            ax_marg_x.hist(r,bins=r_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
        else:
            theta_bins = np.linspace(-90,90,100)
            ax_joint.scatter(angle,centre_mass_z,zorder=zorder,s=1,label=labels[i],marker=markers[i],facecolors="None",edgecolors=colours[i])
            ax_marg_x.hist(angle,bins=theta_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
        ax_marg_y.hist(centre_mass_z,bins=com_bins,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder,density=False)
else:
    facecolours = ['none','red']
    for i,v in enumerate([-1,2]):
        ii=0
        zorder=-i

        if v == -1:
            label = 'Trapped'
        else:
            label = r'$v_{f} =$ ' + str(v)
        #for ii,atom in enumerate(['n','o']):
        #indices = (np.where((np.array(results['final_v'])==v) & (np.array(results['atom_first'])==atom)))[0]
        indices = (np.where((np.array(results['final_v'])==v)))[0]
        
        pos1 = np.array((results['pos1']))[indices,:]
        pos2 = np.array((results['pos2']))[indices,:]

        x1 =  pos1[:,0]
        x2 = pos2[:,0]

        z1 =  pos1[:,1]
        z2 = pos2[:,1]

        y1 = pos1[:,2]
        y2 = pos2[:,2]
        centre_mass_z = (y1*O_mass + y2*N_mass) / (O_mass + N_mass)

        COM = (pos1*O_mass + pos2*N_mass) / (O_mass + N_mass)
        dx = pos1[:,0]-pos2[:,0]
        dy = pos1[:,1]-pos2[:,1]
        dz = pos1[:,2]-pos2[:,2]
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        angle = np.arcsin((y1-y2)/r) * 180/np.pi

        if mode == 'r':
            r_bins = np.linspace(1,1.5,100)
            ax_joint.scatter(r,centre_mass_z,zorder=zorder,s=1,label=label,marker=markers[i],facecolors=facecolours[ii],edgecolors=colours[i])
            ax_marg_x.hist(r,bins=r_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
            ax_marg_y.hist(centre_mass_z,bins=com_bins,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder,density=False)
        elif mode == 'xy':
            xy_bins = np.linspace(-10,10,100)
            ax_joint.scatter(COM[:,0],COM[:,1],zorder=zorder,s=1,label=label,marker=markers[i],facecolors=facecolours[ii],edgecolors=colours[i])
            ax_marg_x.hist(COM[:,0],bins=xy_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
            ax_marg_y.hist(COM[:,1],bins=xy_bins,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder,density=False)
        else:
            theta_bins = np.linspace(-90,90,100)
            ax_joint.scatter(angle,centre_mass_z,zorder=zorder,s=1,label=label,marker=markers[i],facecolors=facecolours[ii],edgecolors=colours[i])
            ax_marg_x.hist(angle,bins=theta_bins,color=colours[i],alpha=0.5,zorder=zorder,density=False)
            ax_marg_y.hist(centre_mass_z,bins=com_bins,color=colours[i],alpha=0.5,orientation="horizontal",zorder=zorder,density=False)
        







#Limits
ax_joint.legend(ncol=2,handletextpad=0.15,columnspacing=0.2,fancybox=True,framealpha=1,handlelength=1,bbox_to_anchor=(0.5, 1.05), loc='center')


annotate_args = {'xy' : (0.01,0.05), 'xycoords' : 'axes fraction'}
ax_joint.annotate(r'$v_i = 3$',ha="left", **annotate_args)
annotate_args['xy'] = (0.01,0.9)
if mode == 'r':
    ax_joint.set_xlim(1.,1.5)
    ax_marg_x.set_xlim(1.,1.5)
    ax_joint.annotate(r'(b)',ha="left", **annotate_args)
    ax_joint.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax_joint.xaxis.set_major_locator(MultipleLocator(0.1))
    ax_joint.set_xlabel(r"r / $\mathrm{\AA{}}$")
    ax_marg_x.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax_marg_x.xaxis.set_major_locator(MultipleLocator(0.2))
    ax_marg_y.set_ylim(1.25,3)
    ax_joint.set_ylim(1.25,3)
    ax_joint.set_ylabel(r"COM height / $\mathrm{\AA{}}$")
elif mode == 'xy':
    ax_joint.set_ylabel(r"y / $\mathrm{\AA{}}$")
    ax_joint.set_xlabel(r"x / $\mathrm{\AA{}}$")
    ax_marg_x.set_xlim(-10,10)
    ax_joint.set_xlim(-10,10)
    ax_marg_y.set_ylim(-1,10)
    ax_joint.set_ylim(-1,10)
else:
    ax_marg_x.set_xlim(-90,90)
    ax_joint.set_xlim(-90,90)
    ax_joint.annotate(r'(a)',ha="left", **annotate_args)
    ax_joint.xaxis.set_minor_locator(MultipleLocator(5))
    ax_joint.xaxis.set_major_locator(MultipleLocator(30))
    ax_joint.set_xlabel(r"$\theta$")
    ax_marg_x.xaxis.set_minor_locator(MultipleLocator(5))
    ax_marg_x.xaxis.set_major_locator(MultipleLocator(30))
    ax_marg_y.set_ylim(1.25,3)
    ax_joint.set_ylim(1.25,3)
    ax_joint.set_ylabel(r"COM height / $\mathrm{\AA{}}$")




ax_joint.yaxis.set_major_locator(MultipleLocator(0.25))
ax_joint.yaxis.set_minor_locator(MultipleLocator(0.05))

ax_marg_y.yaxis.set_major_locator(MultipleLocator(0.25))
ax_marg_y.yaxis.set_minor_locator(MultipleLocator(0.05))

plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
# Set labels on joint
# 



# Set labels on marginals
ax_marg_x.set_ylabel('Frequency')
ax_marg_y.set_xlabel('Frequency')

# 
plt.subplots_adjust(hspace=0.6,wspace=0.6)
fig.set_figwidth(3.25)
fig.set_figheight(3.0)

fig.savefig('scatter.pdf',transparent=True,bbox_inches='tight')
fig.savefig('scatter.png',dpi=300,transparent=True,bbox_inches='tight')