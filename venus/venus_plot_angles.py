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

labels = ['Single bounce', 'Double', 'Multi', 'All bounce']
colours = ['dodgerblue','maroon','navy','black']
linestyles=['-.','--',':','-']
markers =['^','s','.','o']
def cosine_fit(x,m,x0):

    y = np.cos(x-x0)**m

    return y
filenames = sys.argv[2:]
angles = []
bounces = []
ntrajs=0
for i,filename in enumerate(filenames):

    print(filename)
    with open(filename) as f:
        for line in f:

            if 'Trajectory' in line:
                ntrajs += 1
            
            elif 'Lifetime' in line:
                numbers = line.replace(',','')
                theta = float(numbers.split()[6])
                Chi = float(numbers.split()[9])
                
                if Chi > 180 and Chi < 360:
                    scat = -1*theta
                else:
                    scat = theta
                angles.append(scat)
            elif 'bounces' in line:
                bounces.append(int(line.split()[-1]))

angles = np.array(angles) * np.pi/180
bounces = np.array(bounces)

single_bounce_angles = []
double_bounce_angles = []
multi_bounce_angles = []
for t in range(len(bounces)):
            if bounces[t] == 1:
                single_bounce_angles.append(angles[t])

            elif bounces[t] == 2:
                double_bounce_angles.append(angles[t])
            
            elif bounces[t] > 2:
                multi_bounce_angles.append(angles[t])



sb = np.array(single_bounce_angles)
ab = np.array(angles)



# fig, ax = plt.subplots(1, 1, figsize=(3,3.25))
# bins = np.linspace(-110*np.pi/180, 100*np.pi/180, 40)


# x = (bins[1:] + bins[:-1]) / 2

# hist, bin_edges = np.histogram(ab,bins=bins)
# hist = hist / np.sin(x)
# hist = hist / np.max(hist)

# ax.plot(x*180/np.pi,hist,'.-')
# ax.set_xlim(-30,30)
# fig.savefig('angles.pdf',transparent=True,bbox_inches='tight')




fig, ax = plt.subplots(1, 1, figsize=(3,3.25), subplot_kw={'projection': 'polar'})
bins = np.arange(-100,100,10) * np.pi/180


#bin centres
x = (bins[1:] + bins[:-1]) / 2


#sb
hist, bin_edges = np.histogram(sb,bins=bins)
hist = hist / np.abs(np.sin(x))
hist = hist / np.max(hist)
ax.plot(x,hist, '.',color=colours[0],label = labels[0],marker=markers[0],markersize=3,mfc='none')


#ab
hist, bin_edges = np.histogram(ab,bins=bins)
print(bin_edges * 180/np.pi)
hist = hist / np.abs(np.sin(x))
hist = hist / np.max(hist)
ax.plot(x,hist, '.',color=colours[3],label = labels[3],marker=markers[3],markersize=3,mfc='none')



theta = np.linspace(0,2*np.pi,200)
ax.plot(theta,np.cos(theta),linestyle='--',color='firebrick',label=r'cos$(\theta)$',linewidth=.9)

popt, pcov = scipy.optimize.curve_fit(cosine_fit,np.abs(x),hist)
print('m, theta0')
print(popt)
ax.plot(theta,cosine_fit(theta,*popt),'-',color=colours[3],linewidth=0.9,linestyle=linestyles[3],label = r'cos$^{:.0f}(\theta - {:.0f})$'.format(popt[0],popt[1]*180/np.pi))

ax.tick_params(axis='both', which='major')
ax.set_theta_zero_location("N")
ax.set_xticks(np.arange(-90, 100, 30)*np.pi/180)

plt.legend(ncol=2,handletextpad=0.15,columnspacing=0.6,fancybox=True,framealpha=0,handlelength=2,bbox_to_anchor=(0.5, 0.95), loc='center')
ax.set_xlim(-np.pi/2,np.pi/2)
ax.set_ylim(0,1)


plt.annotate("", 
    xy=(0,0.5), xytext=(0,0.9),
    arrowprops=dict(facecolor='black', shrink=0.1, width=0.7, headwidth=6))

ax.set_xlabel('Distribution', labelpad=-20)

ax.set_position( [0.1, -0.45, 0.8, 2.])
fig.savefig('angles.pdf',transparent=True,bbox_inches='tight')
fig.savefig('angles.tiff',transparent=True,bbox_inches='tight',dpi=600)
fig.savefig('angles.eps',transparent=True,bbox_inches='tight')
