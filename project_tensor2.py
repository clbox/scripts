import numpy as np
from sys import argv as argv


#Now reads from friction_tensor.out and NORMALMODES_OUTPUT

lines = open(argv[1]).readlines() #friction_tensor.out
lines2 = open(argv[2]).readlines() #NORMALMODES_OUTPUT


ndim = (len(lines)-1)/2 #dimensions of tensor assuming one title line and one header line for each line
print(ndim)
friction_tensor = np.zeros((ndim,ndim))
i = 0
for line in lines:
    if '#' in line:
        continue
    else:
        for j in range(ndim):
            friction_tensor[i,j]=float(line.split()[j])

        i+=1


modes = np.zeros([ndim,ndim])
for i, line in enumerate(lines2):
  if 'Zero-point energy' in line:
      j = i+1
      for a in range(ndim):
          for b in range(ndim):
              modes[a,b] = \
                      float(lines2[j].split()[b])
          j += 1

for i in range(ndim):
    modes[i,:]/=np.linalg.norm(modes[i,:])

E,V  = np.linalg.eig(friction_tensor[:,:])
E = sorted(E)
string = ' '
for e in E:
    string += str(abs(1./e)) + ' '
string += '\n\n'
for i,e in enumerate(E):
    mode = modes[i,:]/np.linalg.norm(modes[i,:])
    #mode[6:] = 0.
    #mode /= np.linalg.norm(mode)
    string +=  str((1./(np.dot(mode,np.dot(friction_tensor[:,:],mode))))) + ' '
print string

#transform friction matrix into normalmode space
A = np.dot(modes,np.dot(friction_tensor,modes.transpose()))
print 'A'
print A ,'\n'
for i in range(len(E)):
    print 1./A[i,i]
    
fl = open('relax_rate', 'w')
for i in range(len(E)):
    string = '{0:.15f}'.format(A[i,i])
    fl.write(string+'\n')
fl.write('\n')
fl.write('\n')
for i in range(len(E)):
    string = '{0:.15f}'.format(1/A[i,i])
    fl.write(string+'\n')
fl.close()




print 'A trace'
print A.trace(), '\n'
print 'friction_tensor trace'
print friction_tensor.trace(), '\n'

print 'normalmode transformed friction tensor' 
for i in range(ndim):
    string = ''
    for j in range(ndim):
        string += ' {0:14.8f} '.format(A[i,j])
    print string

