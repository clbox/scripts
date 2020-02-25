import numpy as np
from coolvib.tools.spectrum import read_memory_kernel
import matplotlib.pyplot as plt
from ase.units import _hbar, J, s, fs
hbar = _hbar * J
from scipy.interpolate import interp1d


path_to_calc = '/storage/chem/msrvhs/work/prb_16_rerun/co_cu100/c2x2/medium/nvt/calcs/101/friction_memory_kernel.out'
bins,re,im,dimension,max_e = read_memory_kernel(path_to_calc)


frequencies=(bins/hbar)/s

#assume uniformly spaced




vals = np.append(re[0,0,:],re[0,0,1:]*0)
#0 --> positive then negative
frequencies = np.append(frequencies,-frequencies[1:])
idx = np.argsort(frequencies)
times = 1/frequencies

f_interval = frequencies[1]-frequencies[0]
t = np.linspace(0,1/f_interval,len(frequencies))
ft = np.fft.fft(vals)


fig, ax = plt.subplots(6)
ax[0].plot(frequencies[idx],vals[idx],marker='.',markersize=1)
ax[1].plot(t/fs,ft.real,'b-',marker='.')
ax[1].set_xlim(0,5)
ax[2].plot(frequencies[idx],np.fft.fft(ft)[idx])


#Interpolated

new_frequencies = np.linspace(np.min(frequencies),np.max(frequencies),10000)
fit = interp1d(frequencies,vals)
new_vals  = fit(new_frequencies)
idx2 = np.argsort(new_frequencies)
new_f_interval = new_frequencies[1]-new_frequencies[0]
new_t = np.linspace(0,1/new_f_interval,len(new_frequencies))
new_ft = np.fft.fft(new_vals)



ax[3].plot(new_frequencies[idx2],new_vals[idx2],'r-',marker='.',markersize=1)
ax[4].plot(new_t/fs,new_ft.real,'r-',marker='.')
ax[4].set_xlim(0,5)
ax[5].plot(new_frequencies[idx2],np.fft.ifft(new_ft)[idx2],'r-')
plt.show()