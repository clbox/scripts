import numpy as np
import os
import sys


#ARG 1 = friction memory kernel
#ARG 2 = NORMALMODES_OUTPUT (hessian)

#OUTPUT = relaxation_rates_freqs.out

path = sys.argv[1]

#FRICTION MEMORY PARSER
head_count =0
header = ["No of","Discretization","Number of Bins","Excitation energy","==========","k-point","Friction"] #skip lines
with open(path, "r") as f:
    for line in f:
        if "Friction" in line:
            dimension = int(line.split()[3])
            head_count += 1
        if "Discretization" in line:
            discretization=float(line.split()[-1])
        if any(x in line for x in header):
            continue
        max_e = float(line.split()[0])
print("Friction max energy = "+str(max_e))
print("The dimensions of the tensor are " + str(dimension) + "x" + str(dimension))
elements = (((dimension*dimension)-dimension)/2)+dimension
print("There are " + str(elements) + " coupling components")
if elements < head_count:
    n_spin = 2
    print("This system is spin unrestricted")

bins=np.zeros((int(max_e/discretization)+1))
print(len(bins))
re_memory_kernel=np.zeros((dimension,dimension,len(bins)))
im_memory_kernel=np.zeros_like(re_memory_kernel)

with open(path, "r") as f:
    for line in f:
        if "Friction" in line:
            i = int(line.split()[3])
            j = int(line.split()[4])
            head_count += 1
            c=0
        if any(x in line for x in header):
            continue
        else:
            re_memory_kernel[i-1,j-1,c]=float(line.split()[1])
            im_memory_kernel[i-1,j-1,c]=float(line.split()[2])
            bins[c]=float(line.split()[0])
            c +=1
#VIBRATIONAL HESSIAN PARSER
vib_output = sys.argv[2]
lines2 = open(vib_output).readlines()
modes = np.zeros([dimension,dimension])
for i, line in enumerate(lines2):
    if 'Zero-point energy' in line:
        j = i+1 
        for a in range(dimension):
            for b in range(dimension):
                modes[a,b] = float(lines2[j].split()[b])
            j += 1

for i in range(dimension):
    modes[i,:]/=np.linalg.norm(modes[i,:])


#OUTPUT
#TODO output 1xmodes graph as well

output_filename='relaxation_rates_freq.out'

memory_lifetimes_p = np.zeros([dimension,len(bins)])
re_memory_p=np.zeros_like(re_memory_kernel)

for w in range(len(bins)):
    re_memory_p[:,:,w] = np.dot(modes,np.dot((re_memory_kernel[:,:,w]),modes.transpose()))
    for iii in range(dimension):
        memory_lifetimes_p[iii,w] = 1/re_memory_p[iii,iii,w]
    

with open(output_filename,"w") as f:
    f.write('Excitation energy / eV  Vibrational Relaxation rate / ps^-1 \n')
    for iii in range(dimension):
        f.write('Mode = {}\n'.format(iii))
        string = ''
        for w in range(len(bins)):
            string += '{0:.6f}'.format(bins[w])
            string += '      '
            string += '{0:.8f}'.format(re_memory_p[iii,iii,w])
            string += '\n'
        f.write(string)






