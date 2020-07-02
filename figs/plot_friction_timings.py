import numpy as np
import glob
import sys
import os
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.lines import Line2D

filenames = sys.argv[1:]
linewidth = 3
markersize = 12
#annotate=True
annotate=False

fig, ax = plt.subplots(1, 1, sharex='all',sharey='all')#, constrained_layout=True)

original_args = {'marker' : 'o', 'linestyle' : '-',  'alpha' : 1.0, 'mfc' : 'red','linewidth' : linewidth, 'markersize' : markersize }
new_args = {'marker' : 'o','linestyle' : '--', 'alpha' : 1.0,'mfc' :'white', 'linewidth' : linewidth, 'markersize' : markersize }

colours = ['red','blue','purple']

x_scale = [8,16,24,32,40,64]

args = [original_args,new_args]

for i,filename in enumerate(filenames):

    data = np.loadtxt(filename)

    timings = data[:,0]

    scfs = data[:,1]

    n_calcs = int(len(timings) // 3)


    timings = timings.reshape(n_calcs,3)
    scfs = scfs.reshape(n_calcs,3)

    for ii in range(len(colours)):
        (args[0])['mfc'] = colours[ii]

        if annotate and ii<2:
            for c in range(n_calcs):
                ax.annotate(str(int(scfs[c,ii])),xy=(x_scale[c], timings[c,ii]))
        ax.plot(x_scale,timings[:,ii],color=colours[ii],**args[i])

    #(args[0])['mfc'] = 'black'
    #ax.plot(x_scale,np.sum(timings[:,:],axis=1),color='black',**args[i])



ax.set_xticks((x_scale))

ax.set_ylim(0,10000)

font='Times New Roman'
fontsize=22

for tick in ax.get_xticklabels():
    tick.set_fontname(font)
for tick in ax.get_yticklabels():
    tick.set_fontname(font)
ax.tick_params(axis='both', which='major', labelsize=fontsize)

ax.set_ylabel('Total time / s',fontsize=fontsize,fontname=font,color='black')
ax.set_xlabel(r"$N_k$",fontsize=fontsize,fontname=font)

fig.set_figheight(5.0)
fig.set_figwidth(7.0)

plt.gcf().subplots_adjust(left=0.3,bottom=0.2)

lines = []
lines.append(Line2D([0], [0], color='black', linewidth=2,marker='o',mfc='black', linestyle='-'))
lines.append(Line2D([0], [0], color='black', linewidth=2,marker='o',mfc='white', linestyle='--'))

lines2 = [Line2D([0], [0], color=c, linewidth=2, linestyle='-') for c in colours]
for line in lines2:
    lines.append(line)

labels = ['180903','200617','GS SCF', 'FD SCF', r'Calc Tensor']
ax.legend(lines, labels,fancybox=True,framealpha=0,fontsize=12)


fig.savefig('timing.pdf',transparent=True)#,bbox_inches='tight')