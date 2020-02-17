#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import glob
import sys


fig, ax = plt.subplots(2, 1, sharex='all')#, sharey='all')
filename = sys.argv[1]
#column = sys.argv[2]
x = []
y = []

data = np.loadtxt(filename)



#TOTAL ENERGY
y = data[:,4]

ax[0].plot(data[:,0],y, marker='.', markersize=6)
#plt.show()
fig.text(0.01, 0.7, r'Total energy / eV ', va='center', rotation='vertical',fontsize=15)



################################################

#TEMPERATURE

y2 = data[:,1]


ax[1].plot(data[:,0],y2, marker='.', markersize=6)
#plt.show()

fig.text(0.5, 0.07, 'Time / ps', ha='center',fontsize=15)
fig.text(0.01, 0.3, r'Temperature / K ', va='center', rotation='vertical',fontsize=15)


fig.set_figheight(10)
fig.set_figwidth(8)
fig.savefig('nvt_plot.pdf',transparent=True,bbox_inches='tight')
