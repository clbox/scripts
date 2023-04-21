
import numpy as np
import scipy.sparse as sp
import struct
import time


def read_elsi_to_csc(filename):
    mat = open(filename,"rb")
    data = mat.read()
    mat.close()

    # Get header
    start = 0
    end = 64
    header = struct.unpack("i"*16,data[start:end])

    # Number of basis functions (matrix size)
    n_basis = header[3]

    # Total number of non-zero elements
    nnz = header[5]

    # Get column pointer
    start = end
    end = start+n_basis*4
    col_ptr = struct.unpack("i"*n_basis,data[start:end])
    col_ptr += (nnz+1,)
    col_ptr = np.array(col_ptr)

    # Get row index
    start = end
    end = start+nnz*4
    row_idx = struct.unpack("i"*nnz,data[start:end])
    row_idx = np.array(row_idx)

    # Get non-zero value
    start = end

    if header[2] == 0:
        # Real case
        end = start+nnz*8
        nnz_val = struct.unpack("d"*nnz,data[start:end])
    else:
        # Complex case
        end = start+nnz*16
        nnz_val = struct.unpack("d"*nnz*2,data[start:end])
        nnz_val_real = np.array(nnz_val[0::2])
        nnz_val_imag = np.array(nnz_val[1::2])
        nnz_val = nnz_val_real + 1j*nnz_val_imag

    nnz_val = np.array(nnz_val)

    # Change convention
    for i_val in range(nnz):
        row_idx[i_val] -= 1

    for i_col in range(n_basis+1):
        col_ptr[i_col] -= 1

    return sp.csc_matrix((nnz_val,row_idx,col_ptr),shape=(n_basis,n_basis))


def parse_chem_pot(aims_file):
    chem_pot = 0
    with open(aims_file, "r") as af:
        for line in af:
            if '**FRICTION**' in line:
                break
            if '| Chemical potential (Fermi level):' in line:
                chem_pot = float(line.split()[-2])
    return chem_pot # eV

def parse_evs(evs_file):
    with open(evs_file,"r") as f:
        for i in range(3):
            line = f.readline()
        line = f.readline()
        ncols = len(line.split())


    evs = np.loadtxt(evs_file,skiprows=3,usecols=range(1,ncols))
    # print(np.shape(evs))

    return evs


def find_evs_sensible_bounds(evs,chem_pot,max_energy_from_fermi):
    # We find indices for upper and lower bounds of eigenvalues
    # to reduce size of problem
    min_bounds = np.zeros((np.shape(evs)[0]),dtype=int)
    max_bounds = np.zeros((np.shape(evs)[0]),dtype=int)


    for i_k_point in range(np.shape(evs)[0]):
        for i,ev in enumerate(evs[i_k_point,:]):
            if ev < (chem_pot-max_energy_from_fermi):
                min_bounds[i_k_point]= i
                continue
            elif ev > (chem_pot+max_energy_from_fermi):
                max_bounds[i_k_point]= i
                break
    return min_bounds,max_bounds

def gaussian_function(x,x0,s):
    return 1./np.sqrt(2*np.pi*s*s)*np.exp(-0.5*(((x-x0) / (s))**2))


def parse_timestep_data_kpoint(dirname,i_k_point,i_atom,i_cart,i_spin,parse_evecs=True):


    aims_filename = dirname+"aims.out"
    evs_filename = dirname+"friction_KS_eigenvalues.out"
    
    ham1_filename = dirname+"first_order_H_atom_ATOMID_cart_CARTID_k_KID.csc"
    ovlp1_filename = dirname+"first_order_S_atom_ATOMID_cart_CARTID_k_KID.csc"

    # Parse Fermi level (chemical potential)
    # start = time.time() 
    chem_pot = parse_chem_pot(aims_filename)
    # end = time.time()
    # print('        Time for 1 parse chempot/ s: '+str(end - start))
    # print('Chemical potential: '+str(chem_pot)+' / eV')

    evs_file = evs_filename
  
    ham1_file = ham1_filename.replace('ATOMID','{:06d}'.format(i_atom+1)).replace('CARTID','{:1d}'.format(i_cart+1)).replace('KID','{:06d}'.format(i_k_point+1))
    ovlp1_file = ovlp1_filename.replace('ATOMID','{:06d}'.format(i_atom+1)).replace('CARTID','{:1d}'.format(i_cart+1)).replace('KID','{:06d}'.format(i_k_point+1))

    # start = time.time() 
    evs = parse_evs(evs_file)
    # end = time.time()
    # print('        Time for 1 parse ev/ s: '+str(end - start))

    # start = time.time()
    ham1= read_elsi_to_csc(ham1_file)
    ovlp1 = read_elsi_to_csc(ovlp1_file)
    # end = time.time()
    # print('        Time for 1 parse ham1,ovlp1/ s: '+str(end - start))


    if parse_evecs:
        evecs_filename = dirname+"C_spin_SPINID_kpt_KID.csc"
        evecs_file = evecs_filename.replace('SPINID','{:02d}'.format(i_spin+1)).replace('KID','{:06d}'.format(i_k_point+1))
        evecs = read_elsi_to_csc(evecs_file)
    else:
        evecs = 0.

    return chem_pot, evs, evecs, ham1, ovlp1


def parse_epc_data_kpoint(dirname,i_k_point,i_atom,i_cart,i_spin):


    # aims_filename = dirname+"aims.out"
    # evs_filename = dirname+"friction_KS_eigenvalues.out"
    
    epc_filename = dirname+"epc_atom_ATOMID_cart_CARTID_k_KID.csc"


    # Parse Fermi level (chemical potential)
    # start = time.time() 
    # chem_pot = parse_chem_pot(aims_filename)
    # end = time.time()
    # print('        Time for 1 parse chempot/ s: '+str(end - start))
    # print('Chemical potential: '+str(chem_pot)+' / eV')

    # evs_file = evs_filename
  
    epc_file = epc_filename.replace('ATOMID','{:06d}'.format(i_atom+1)).replace('CARTID','{:1d}'.format(i_cart+1)).replace('KID','{:06d}'.format(i_k_point+1))


    # start = time.time() 
    # evs = parse_evs(evs_file)
    # end = time.time()
    # print('        Time for 1 parse ev/ s: '+str(end - start))

    # start = time.time()
    epc= read_elsi_to_csc(epc_file)

    # end = time.time()
    # print('        Time for 1 parse epc/ s: '+str(end - start))


   

    return epc 

def parse_ev_data(dirname):


    aims_filename = dirname+"aims.out"
    chem_pot = parse_chem_pot(aims_filename)

    evs_filename = dirname+"friction_KS_eigenvalues.out"
    evs_file = evs_filename
    evs = parse_evs(evs_file)


    return chem_pot, evs