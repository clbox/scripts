import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                       AutoMinorLocator)

fig4, ax = plt.subplots(1, 1, sharex='all')#, sharey='all')

color_idx = np.linspace(0, 1, len(sys.argv[1:]))
max_val = 0
min_val = 0

for t, (filename, i) in enumerate(zip(sys.argv[1:], color_idx)):
    print(filename)
    a = np.loadtxt(filename,skiprows=4)
    ax.plot(a[:,0],a[:,1],'-',linewidth = 1.5,label=str(filename),
            color=plt.cm.cool(float(t)/float(len(sys.argv[1:]))))
    ax.set_xlim(0,1)
    index = int(np.where(a[:,0]==1)[0])
    if np.max((a[:,1])[:index]) > max_val:
        max_val = np.max((a[:,1])[:index])
    if np.min((a[:,1])[:index]) < min_val:
        min_value = np.min((a[:,1])[:index])
    ax.set_ylim(min_val,max_val)

ax.legend()
fig4.set_figheight(6)
fig4.set_figwidth(6)
fig4.text(0.5, 0.01, "Excitation energy / eV", ha='center',fontsize=15)
fig4.text(0.01, 0.5, r'JDOS / eV$^{-2}$', va='center', rotation='vertical',fontsize=15)
fig4.savefig('jdos.pdf',transparent=True,bbox_inches='tight')
