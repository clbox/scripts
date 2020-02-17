import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import glob
import sys


filename = sys.argv[1]
#column = sys.argv[2]
x = []
y = []

data = np.loadtxt(filename)

plt.figure(figsize=(3,3))


y = data[:,4]

plt.plot(data[:,0],y, marker='.', markersize=6)
#plt.show()
plt.xlabel('Time / ps')
plt.ylabel('Energy / eV')

#plt.savefig(str(filename)+'_E.pgf',bbox_inches='tight')
plt.savefig(str(filename)+'_E.pdf',bbox_inches='tight')



################################################


plt.figure(figsize=(3,3))

y2 = data[:,1]


plt.plot(data[:,0],y2, marker='.', markersize=6)
#plt.show()
plt.xlabel('Time / ps')
plt.ylabel('Temperature / T')

#plt.savefig(str(filename)+'_T.pgf',bbox_inches='tight')
plt.savefig(str(filename)+'_T.pdf',bbox_inches='tight')

