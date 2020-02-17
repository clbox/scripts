#Reads in tensor from aims.out. Made specifically for large tensors where each row of the tensor is
#printed over multiple lines in aims.out but should work in normal (small) tensors.

import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import sys

#Reading tensor from aims.out
lines = open(sys.argv[1]).readlines()
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
friction_tensor = np.zeros((ndim,ndim))
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


b= (np.ndarray.flatten(friction_tensor))
print(' '.join(map(str, b))) 
