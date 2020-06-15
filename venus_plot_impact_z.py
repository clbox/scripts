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

filenames = sys.argv[1:]


labels = ['Trapped',r'$\nu_f = 1$',r'$\nu_f = 2$']

colours = ['blue','red','orange']
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