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


filename = sys.argv[1]

outfile = filename.replace('.txt','')
# Read the array from disk
impact_geos = np.loadtxt(filename)

# Note that this returned a 2D array!
print(impact_geos.shape)

impact_geos = impact_geos.reshape(impact_geos.shape[0]//2,2,3)

print(impact_geos.shape)


fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')


ax.scatter(impact_geos[:,0,0],impact_geos[:,0,2],c='red',s=20,alpha=0.7,label='Oxygen')
ax.scatter(impact_geos[:,1,0],impact_geos[:,1,2],c='blue',s=20,alpha=0.7,label='Nitrogen')

ax.set_xlim(-5,7)
ax.set_ylim(0,7)
lgnd = ax.legend(loc=1)
lgnd.legendHandles[0]._sizes = [200]
lgnd.legendHandles[1]._sizes = [200]

fig.set_figheight(4)
fig.set_figwidth(5)
fig.text(0.5, 0.00, r"$x$ / $\AA{}$", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'$z$ / $\AA{}$', va='center', rotation='vertical',fontsize=15)
fig.savefig(outfile+'.pdf',transparent=True,bbox_inches='tight')



############ With histograms ################ 
x1 = impact_geos[:,0,0] #Ox
y1 = impact_geos[:,0,2] #Oz

x2 = impact_geos[:,1,0] #Nx
y2 = impact_geos[:,1,2] #Nz

fig = plt.figure()

gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])

bins_x = np.linspace(-5, 7, 100)
bins_y = np.linspace(0, 7, 100)

ax_joint.scatter(x1,y1,c='red',s=20,alpha=0.7,label='Oxygen')
ax_joint.scatter(x2,y2,c='blue',s=20,alpha=0.7,label='Nitrogen')
ax_marg_x.hist(x1,bins_x,color='red',alpha=0.5)
ax_marg_x.hist(x2,bins_x,color='blue',alpha=0.5)
ax_marg_y.hist(y1,bins_y,color='red',alpha=0.5,orientation="horizontal")
ax_marg_y.hist(y2,bins_y,color='blue',alpha=0.5,orientation="horizontal")

#Limits
ax_marg_y.set_ylim(0.5,4.5)
ax_marg_x.set_xlim(-5,7)
ax_joint.set_xlim(-5,7)
ax_joint.set_ylim(0.5,4.5)

ax_joint.legend(loc=1)

# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel(r'$x$ / $\AA{}$')
ax_joint.set_ylabel(r'$z$ / $\AA{}$')

# Set labels on marginals
ax_marg_y.set_xlabel('Frequency')
ax_marg_x.set_ylabel('Frequency')

fig.savefig(outfile+'_hist.pdf',transparent=True,bbox_inches='tight')
fig.savefig(outfile+'_hist.png',dpi=300,transparent=True,bbox_inches='tight')