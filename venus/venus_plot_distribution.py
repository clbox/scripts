import numpy as np
from ase.units import _hbar, J, s, fs
from ase import Atoms
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, MaxNLocator)
import sys
import scipy
from scipy import optimize


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

matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"


#mode 0 = bomd, mode 1 = ldfa, mode 2 = tdpt
mode = int(sys.argv[1])

filenames = sys.argv[2:]

nstates = len(filenames)

#PLOT TDPT ENERGY LOSS 
if mode == 2:
    fig, ax = plt.subplots(nstates, 1, sharex='all')#,sharey='all')

#PLOT ENERGY REDISTRIBUTION

ymax=2.5

empty_ones = []
ntraj_list = []
ntraj_single_list = []
ntraj_double_list = []
ntraj_multi_list = []
state_list=[]
e_diff_list=[]
angles = []
single_bounce_angles = []
double_bounce_angles = []
multi_bounce_angles = []

orientated=False
n_first_list = []
o_first_list = []


n_first_single_list = []
o_first_single_list = []
for i,filename in enumerate(filenames):

    print(filename)

    ntrajs = 0
    ntraj_single = 0
    ntraj_double = 0
    ntraj_multi = 0
    ntraj_nfirst = 0
    ntraj_ofirst = 0
    ntraj_nfirst_single = 0
    ntraj_ofirst_single = 0
    misc = []
    abs_e = []
    per_e = []
    init_e = []
    final_e = []
    n_bounces = []
    atom_first = []

    with open(filename) as f:
        for line in f:

            if 'Trajectory' in line:
                ntrajs += 1
            
            elif 'Lifetime' in line:
                numbers = line.replace(',','')
                lifetime = float(numbers.split()[2])
                scat = float(numbers.split()[7])
                jf = int(numbers.split()[-1])

                misc.append([lifetime,scat,jf])

            elif 'Total' in line:
                numbers = line.replace(',','')
                d = float(numbers.split()[7])
                theta = float(numbers.split()[10])
                phi = float(numbers.split()[13])
                X = float(numbers.split()[16])
                Y = float(numbers.split()[19])
                Z = float(numbers.split()[22])
                tot = float(numbers.split()[-1])

                abs_e.append([d,theta,phi,X,Y,Z,tot])

            elif '%' in line:
                numbers = line.replace(',','')
                d = float(numbers.split()[6])
                theta = float(numbers.split()[9])
                phi = float(numbers.split()[12])
                X = float(numbers.split()[15])
                Y = float(numbers.split()[18])
                Z = float(numbers.split()[21])

                per_e.append([d,theta,phi,X,Y,Z])

            elif 'Initial' in line:
                numbers = line.replace(',',' ')
                i_v = float(numbers.split()[8])
                i_r = float(numbers.split()[9])
                i_t = float(numbers.split()[10])

                init_e.append([i_v,i_r,i_t])

            elif 'Final' in line:
                numbers = line.replace(',',' ')
                f_v = float(numbers.split()[8])
                f_r = float(numbers.split()[9])
                f_t = float(numbers.split()[10])

                final_e.append([f_v,f_r,f_t])

            elif 'bounces' in line:
                n_bounces.append(int(line.split()[-1]))

            elif 'first' in line:
                orientated=True
                if 'N first' in  line:
                    ntraj_nfirst+=1
                if 'O first' in  line:
                    ntraj_ofirst+=1

                atom_first.append(line.split()[0])


    if 'trapped' not in filename:
        for a,atom in enumerate(atom_first):
            if n_bounces[a] == 1:
                if atom =='N':
                    ntraj_nfirst_single+=1
                if atom =='O':
                    ntraj_ofirst_single+=1


    misc = np.array(misc)
    abs_e = np.array(abs_e)
    per_e = np.array(per_e)
    init_e = np.array(init_e)
    final_e = np.array(final_e)


    if 'trapped' in filename:
        if mode == 2:
            ax[0].boxplot(abs_e,showfliers=False)
            ax[0].text(s=filename,x=2,y=np.max(abs_e)-0.5*np.max(abs_e))
            ax[0].text(s='ntrajs = '+str(ntrajs),x=5,y=np.max(abs_e)-0.5*np.max(abs_e))
            ax[0].set_ylim(0,np.max(abs_e))

        n_trapped = ntrajs

    else:
        if mode == 2:
            ax[i+1].set_ylim(0,ymax)
            ax[i+1].boxplot(abs_e,showfliers=False)
            ax[i+1].text(s=filename,x=2,y=ymax-0.5*ymax)
            ax[i+1].text(s='ntrajs = '+str(ntrajs),x=5,y=ymax-0.5*ymax)

        
        e_diff = final_e - init_e
        print('Average lifetime / fs ' +str(np.average(misc[:,0])))
        print('Average scattering angle ' + str(np.average(misc[:,1])))
        print('Average final rotational state ' + str(np.average(misc[:,2])))
        ntraj_list.append(ntrajs)
        state_list.append(int(filename.split('.')[0]))
        e_diff_list.append(e_diff)

        for t in range(len(n_bounces)):
            if n_bounces[t] == 1:
                ntraj_single += 1
                single_bounce_angles.append(misc[t,1])

            elif n_bounces[t] == 2:
                ntraj_double += 1
                double_bounce_angles.append(misc[t,1])
            
            elif n_bounces[t] > 2:
                ntraj_multi += 1
                multi_bounce_angles.append(misc[t,1])

        ntraj_single_list.append(ntraj_single)
        ntraj_double_list.append(ntraj_double)
        ntraj_multi_list.append(ntraj_multi)

        n_first_list.append(ntraj_nfirst)
        o_first_list.append(ntraj_ofirst)

        n_first_single_list.append(ntraj_nfirst_single)
        o_first_single_list.append(ntraj_ofirst_single)
        angles.append(misc[:,1])

if mode == 2:
    plt.xticks([1,2,3,4,5,6,7],['d',r'$\theta$',r'$\phi$','X','Y','Z','Total'])
    fig.set_figheight(10*nstates/10)
    fig.set_figwidth(7)
    fig.text(0.5, 0.05, "Mode", ha='center',fontsize=15)
    fig.text(0.01, 0.5, 'Energy loss / eV', va='center', rotation='vertical',fontsize=15)
    fig.savefig('summary.pdf',transparent=True,bbox_inches='tight')


print(nstates)
fig2, ax2 = plt.subplots(nstates-1, 1, sharex='all',sharey='all')
for i,filename in enumerate(filenames):
    if 'trapped' in filename:
        continue
    else:
        if nstates-1 == 1:
            ax2.set_xlim(0.5,3.5)
            # ax2[i].set_ylim(0,ymax)
            ax2.boxplot(e_diff_list[i],showfliers=False)
            ax2.text(s=filename,x=0.5,y=0.8)
            ax2.text(s='ntrajs = '+str(ntraj_list[i]),x=0.5,y=0.5)
        else:
            ax2[i].set_xlim(0.5,3.5)
            # ax2[i].set_ylim(0,ymax)
            ax2[i].boxplot(e_diff_list[i],showfliers=False)
            ax2[i].text(s=filename,x=0.5,y=0.8)
            ax2[i].text(s='ntrajs = '+str(ntraj_list[i]),x=0.5,y=0.5)

#Write energy distributions to file
e_diff_list = np.array(e_diff_list)
np.save('e_diff_list.npy',e_diff_list)



plt.xticks([1,2,3],['Vibrational','Rotational','Translational'])
fig2.set_figheight(10*nstates/10)
fig2.set_figwidth(7)
fig2.text(0.5, 0.05, "Component", ha='center',fontsize=15)
fig2.text(0.01, 0.5, r'E$_i$ - E$_f$ / eV', va='center', rotation='vertical',fontsize=15)
fig2.savefig('energy_distribution.pdf',transparent=True,bbox_inches='tight')

ntraj_list = np.array(ntraj_list)

np.savetxt('states.txt', np.c_[state_list,ntraj_list/np.sum(ntraj_list)],fmt='%1.7f')


if not mode==-1:
    np.savetxt('states_1.txt', np.c_[state_list,ntraj_single_list/np.sum(ntraj_single_list)],fmt='%1.7f')
    np.savetxt('states_2.txt', np.c_[state_list,ntraj_double_list/np.sum(ntraj_double_list)],fmt='%1.7f')
    np.savetxt('states_multi.txt', np.c_[state_list,ntraj_multi_list/np.sum(ntraj_multi_list)],fmt='%1.7f')

    if orientated==True:
        np.savetxt('states_n.txt', np.c_[state_list,n_first_list/np.sum(n_first_list)],fmt='%1.7f')
        np.savetxt('states_o.txt', np.c_[state_list,o_first_list/np.sum(o_first_list)],fmt='%1.7f')


        np.savetxt('states_1_n.txt', np.c_[state_list,n_first_single_list/np.sum(n_first_single_list)],fmt='%1.7f')
        np.savetxt('states_1_o.txt', np.c_[state_list,o_first_single_list/np.sum(o_first_single_list)],fmt='%1.7f')


    bounce_prob = []
    bounce_prob.append(np.sum(ntraj_single_list))
    bounce_prob.append(np.sum(ntraj_double_list))
    bounce_prob.append(np.sum(ntraj_multi_list))

    np.savetxt('absolute_bounce.txt',bounce_prob)

state_list.append(-1)
ntraj_list = np.append(ntraj_list,n_trapped)

np.savetxt('absolute_pops.txt',np.c_[state_list,ntraj_list],fmt='%1.6f')

#PLOT SCATTERING ANGULAR DISTRIBUTION
#fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='polar')
fig, ax = plt.subplots(1, 1, figsize=(3,3.25), subplot_kw={'projection': 'polar'})
def cosine_fit(x,m,x0):

    y = np.cos(x-x0)**m

    return y


bins = np.linspace(0, 100, 20)
theta = np.linspace(0,2*np.pi,200)
rho = np.cos(theta)
x = (bins[1:] + bins[:-1]) / 2
labels = ['Single bounce', 'Double', 'Multi', 'All bounce']
colours = ['dodgerblue','maroon','navy','black']
linestyles=['-.','--',':','-']
all_angles = single_bounce_angles+double_bounce_angles+multi_bounce_angles
print(np.shape(all_angles))
list_of_angles = [single_bounce_angles,double_bounce_angles,multi_bounce_angles,all_angles]

markers =['^','s','.','o']
for i in range(4):
    #digitized = np.digitize(list_of_angles[i], bins)
    if i in [1,2]:
        continue
    hist, bin_edges = np.histogram(list_of_angles[i],bins=bins)
    hist = hist / np.max(hist)



    ax.plot(x*np.pi/180,hist, '.',color=colours[i],label = labels[i],marker=markers[i],markersize=3,mfc='none')
    
    popt, pcov = scipy.optimize.curve_fit(cosine_fit,x*np.pi/180,hist)
    print('m, theta0')
    print(popt)
    if i == 3:
        ax.plot(theta,cosine_fit(theta,*popt),'-',color=colours[i],linewidth=0.9,linestyle=linestyles[i],label = r'cos$^{:.0f}(\theta - {:.0f})$'.format(popt[0],popt[1]*180/np.pi))
ax.plot(theta,rho,linestyle='--',color='firebrick',label=r'cos$(\theta)$',linewidth=.9)

ax.tick_params(axis='both', which='major')
ax.set_theta_zero_location("N")
ax.set_xticks(np.arange(-90, 100, 10)*np.pi/180)
#ax.set_yticks([])
#ax.xaxis.set_minor_locator(MultipleLocator(5*np.pi/180))
#ax.legend(loc=8,ncol=3,fancybox=True,framealpha=0)
plt.legend(ncol=2,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 0.95), loc='center')
ax.set_xlim(-np.pi/2,np.pi/2)
ax.set_ylim(0,1)

#plt.gcf().subplots_adjust(bottom=-5.5,top=-1.6)
# fig.set_figheight(2.3)
# fig.set_figwidth(3.25)
ax.set_xlabel('Distribution', labelpad=-20)
# ax.set_xlabel("Scattering angle",fontsize=12,fontname=font)
# ax.set_ylabel(r'Density',fontsize=12,fontname=font)
ax.set_position( [0.1, -0.45, 0.8, 2.])
fig.savefig('angles.pdf',transparent=True,bbox_inches='tight')
fig.savefig('angles.tiff',transparent=True,bbox_inches='tight',dpi=600)
fig.savefig('angles.eps',transparent=True,bbox_inches='tight')
