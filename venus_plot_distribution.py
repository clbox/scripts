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
x16_exp = np.arange(0,17,1)
v16_exp = [0.0,0.0,0.04,0.08,0.13,0.15,0.19,0.11,0.12,0.07,0.04,0.02,0.03,0.02,0.01,0.02,0.02]

x3_exp = np.arange(0,7,1)
v3_exp = [0.0,0.2195652173913044,0.42391304347826086,0.35434782608695653,0.0,0,0]

x15_exp = np.arange(5,16)
v15_exp = [0.115,0.1339,0.194,0.192,0.125,0.082,0.04,0.05,0.019,0.015,0.036]

x11_exp = np.arange(2,12)
v11_exp = [0.047,0.099,0.18,0.18,0.18,0.12,0.068,0.052,0.036,0.025]

x2_exp = np.arange(1,4,1)
v2_exp = [0.33,0.66,0.00038]


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
for i,filename in enumerate(filenames):

    print(filename)

    ntrajs = 0
    ntraj_single = 0
    ntraj_double = 0
    ntraj_multi = 0
    misc = []
    abs_e = []
    per_e = []
    init_e = []
    final_e = []
    n_bounces = []

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
fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
ax.plot(state_list,ntraj_list/np.sum(ntraj_list),
   # '.-',color='purple',label=r'TDPT')# $\times 2$')
    marker='^',linestyle='-',color='red',label=r'BOMD')
    # marker='s',linestyle='-',color='blue',label=r'LDFA $\times 4$')


np.savetxt('states.txt', np.c_[state_list,ntraj_list/np.sum(ntraj_list)],fmt='%1.3f')


if not mode==-1:
    np.savetxt('states_1.txt', np.c_[state_list,ntraj_single_list/np.sum(ntraj_single_list)],fmt='%1.3f')
    np.savetxt('states_2.txt', np.c_[state_list,ntraj_double_list/np.sum(ntraj_double_list)],fmt='%1.3f')
    np.savetxt('states_multi.txt', np.c_[state_list,ntraj_multi_list/np.sum(ntraj_multi_list)],fmt='%1.3f')


    bounce_prob = []
    bounce_prob.append(np.sum(ntraj_single_list))
    bounce_prob.append(np.sum(ntraj_double_list))
    bounce_prob.append(np.sum(ntraj_multi_list))

    np.savetxt('absolute_bounce.txt',bounce_prob)

state_list.append(-1)
ntraj_list = np.append(ntraj_list,n_trapped)

np.savetxt('absolute_pops.txt',np.c_[state_list,ntraj_list],fmt='%1.3f')
#V16
# ax.bar(x16_exp,v16_exp,color='black',label=r'$\nu_i=16$ exp')
# ax.set_ylim(0,0.5)
# ax.set_xlim(0,20)

#V3
# ax.bar(x3_exp,v3_exp,color='black',label=r'$\nu_i=3$ exp')
# ax.set_xlim(0,6)
# ax.set_ylim(0,0.8)

#V15
# ax.bar(x15_exp,v15_exp,color='black',label=r'$\nu_i=15$ exp')
# ax.set_xlim(0,17)
# ax.set_ylim(0,0.8)

#V11
# ax.bar(x11_exp,v11_exp,color='black',label=r'$\nu_i=11$ exp')
# ax.set_xlim(0,20)
# ax.set_ylim(0,0.5)

#V2
ax.bar(x2_exp,v2_exp,color='black',label=r'$\nu_i=2$ exp')
ax.set_xlim(0,4)
ax.set_ylim(0,0.8)


ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.legend()


fig.set_figheight(4)
fig.set_figwidth(5)
fig.text(0.5, 0.00, r"Final vibrational state ($\nu_f$)", ha='center',fontsize=15)
fig.text(0.01, 0.5, 'Population', va='center', rotation='vertical',fontsize=15)
fig.savefig('probability.pdf',transparent=True,bbox_inches='tight')





#PLOT SCATTERING ANGULAR DISTRIBUTION
#fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')



def cosine_fit(x,m,x0):

    y = np.cos(x-x0)**m

    return y


bins = np.linspace(0, 100, 20)
theta = np.linspace(0,2*np.pi,200)
rho = np.cos(theta)
x = (bins[1:] + bins[:-1]) / 2
# for i,state in enumerate(angles):
#         if i > 3:
#             continue
#         sns.distplot(state[:], bins = bins, hist = False, kde = True,
#                  kde_kws = {'shade': False, 'linewidth': 1}, 
#                   label = str(i))

labels = ['Single', 'Double', 'Multi', 'Total']
colours = ['dodgerblue','maroon','navy','black']
all_angles = single_bounce_angles+double_bounce_angles+multi_bounce_angles
print(np.shape(all_angles))
list_of_angles = [single_bounce_angles,double_bounce_angles,multi_bounce_angles,all_angles]


for i in range(4):
    #digitized = np.digitize(list_of_angles[i], bins)
    if i in [1,2]:
        continue
    hist, bin_edges = np.histogram(list_of_angles[i],bins=bins)
    hist = hist / np.max(hist)



    ax.plot(x*np.pi/180,hist, '.',color=colours[i],label = labels[i])
    
    popt, pcov = scipy.optimize.curve_fit(cosine_fit,x*np.pi/180,hist)
    print(popt)
    ax.plot(theta,cosine_fit(theta,*popt),'-',color=colours[i])

    #ax.text(0,1.3+(0.05*i),r'$m = {:.2f}, \theta_0 = {:.2f}$'.format(popt[0],popt[1]*180/np.pi),color = colours[i])
    #sns.distplot(list_of_angles[i], bins = bins, hist = False, kde = True,
    #              kde_kws = {'shade': False, 'linewidth': 1}, 
    #               label = labels[i])
# plt.polar(bins/np.pi,np.cos(bins*np.pi/180)/100,'black')
# plt.polar(bins/np.pi,(np.cos(bins*np.pi/180)**14)/100,'black',linestyle=':')

#plt.polar(theta,rho)
ax.plot(theta,rho,linestyle='--',color='black',label=r'$cos(\theta)$')
font='Arial'
for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)

ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_theta_zero_location("N")
ax.set_xticks(np.arange(-90, 100, 10)*np.pi/180)
ax.set_yticks([])
#ax.xaxis.set_minor_locator(MultipleLocator(5*np.pi/180))
ax.legend(loc=8,ncol=3,fontsize=12,fancybox=True,framealpha=0)
ax.set_xlim(-np.pi/2,np.pi/2)
ax.set_ylim(0,1)
plt.gcf().subplots_adjust(bottom=-1)
fig.set_figheight(2.7)
fig.set_figwidth(3.25)
# ax.set_xlabel("Scattering angle",fontsize=12,fontname=font)
# ax.set_ylabel(r'Density',fontsize=12,fontname=font)
fig.savefig('angles.pdf',transparent=True,bbox_inches='tight')
