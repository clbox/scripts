#Reads in tensor from aims.out. Made specifically for large tensors where each row of the tensor is
#printed over multiple lines in aims.out but should work in normal (small) tensors. With a geometry.in file present
#it labels the tensor with the atoms and the corresponding .csv file. No arguments needed

import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import seaborn as sns; sns.set()

#Reading tensor from aims.out
lines = open('aims.out').readlines()
raw_lines = []
copy = False
for line in lines:
    if '********Printing Friction Tensor in 1/ps*********' in line:
        copy = True
        continue
    if '**********END Printing Friction Tensor************' in line:
        copy = False
    if copy:
        raw_lines.append(line)
        #for j in range(ndim):
            #friction_tensor[i,j]=float(line.split()[j])

lines = open('geometry.in').readlines()

friction_indices = []
for ii,line in enumerate(lines):
    if 'calculate_friction .true.' in line:
        friction_indices.append(ii-1)

friction_atoms = []
for ii,line in enumerate(lines):
    if ii in friction_indices:
        friction_atoms.append(line.split()[-1])
a = len(raw_lines[0].split())
lines_per_line = 0
for i,line in enumerate(raw_lines):
    if i == 0:
        continue
    if len(line.split()) != a:
        lines_per_line += 1
    if len(line.split()) == a:
        break

ndim = (a - 1)*(lines_per_line+1)
print(ndim)
friction_tensor = np.zeros((ndim,ndim))
print('lines per line '+str(lines_per_line+1))
c=0
i=-1
for ii, line in enumerate(raw_lines):
    
    if c == lines_per_line+1 or c==0:
        i+=1
        c=0
    
    for j in range(ndim/(lines_per_line+1)):
        if c ==0:
            friction_tensor[i,j]=float(line.split()[j+1])
        else:
            friction_tensor[i,j+((ndim/(lines_per_line+1))*c)]=float(line.split()[j])
    c+=1

labels = []
for atom in friction_atoms:    
    for cart in ['x','y','z']:
        labels.append(r'$\mathrm{'+atom+'_'+cart+'}$')
ft = pd.DataFrame(friction_tensor,columns=labels,index=labels)
ax = sns.heatmap(ft,linewidth=0,vmin=-0.25,vmax=0.4,center=0,cmap="RdBu",xticklabels=True, yticklabels=True)
ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 10)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, fontsize = 10)
#ax.xaxis.set_ticks_position('top')
ax.get_xaxis().set_label_coords(0.5,1.1)
fig = ax.get_figure()
fig.savefig('tensor_mat.pdf') 
