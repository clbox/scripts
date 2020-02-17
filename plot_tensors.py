#Reads in tensor from aims.out. Made specifically for large tensors where each row of the tensor is
#printed over multiple lines in aims.out but should work in normal (small) tensors. With a geometry.in file present
#it labels the tensor with the atoms and the corresponding .csv file. No arguments needed

import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sys

tensors = [None]*len(sys.argv[1:])
dirs = [None]*len(sys.argv[1:])
for idx,output in enumerate(sys.argv[1:]):
    print(idx)
    dirs[idx]=float((output.split('/'))[0])
    #Reading tensor from aims.out
    lines = open(str(output)).readlines()
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
    
    tensors[idx]=friction_tensor


tensors_ary = np.asarray(tensors)
fig, ax = plt.subplots(ndim, ndim, sharex='all', sharey='all')

for i in range(ndim):
    for j in range(ndim):
        ax[i,j].plot(dirs,tensors_ary[:,i,j],'s-')
        ax[i,j].set_ylim(1.5,1.8)

fig.set_figheight(10)
fig.set_figwidth(10)
fig.text(0.5, 0.01, "Temperature / K", ha='center',fontsize=15)
fig.text(0.01, 0.5, r'$\Lambda\ /\ \mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=15)
fig.savefig('tensor_plot.pdf',transparent=True,bbox_inches='tight')
