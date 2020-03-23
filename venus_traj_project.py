import numpy as np
from ase.units import _hbar, J, s, fs
from ase import Atoms
import glob
from pathlib import Path
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
hbar = _hbar * J * s 
ps = fs*1000
labels = [r'$d$',r'$\phi$',r'$\theta$',r'$X$',r'$Y$',r'$Z$']
line_settings = {
    #'marker' : 's',
    'marker' : None,
    'alpha' : 0.8,
    'markersize' : 2
    }

def get_natoms(traj_no):

    filename = 'fort.{}'.format(1000+traj_no)

    with open(filename) as f:
        n_atoms = int(f.readline().split()[0])
    
    return n_atoms

def parse_traj_summary(traj_no):
    """read fort.26 and gets all relevant trajectory information"""

    filename = 'fort.26'

    with open(filename) as f:
        for line in f:
            if 'TRAJ:' in line:
                if line.split()[1] == str(traj_no):
                    if 'LIFETIME' in line:
                        #fs
                        lifetime = float(line.split()[-2])
                    
                    elif 'PRODUCT IS A DIATOM' in line:
                        Nf = float(line.split()[-3]) #vib
                        Jf = float(line.split()[-1]) #rot

                    elif 'SCATTERING_ANG' in line:
                        scat_angle = float(line.split()[-2])
                        break
                elif int(line.split()[1]) > traj_no:
                    break


    return lifetime*fs,Nf,Jf,scat_angle

def parse_traj_velocities(traj_no):
    """reads in fort.(1000+traj_no) and returns velocities for all atoms and 
    time steps"""

    filename = 'fort.{}'.format(1000+traj_no)

    with open(filename) as f:
        n_atoms = get_natoms(traj_no)
        lines_per_traj = n_atoms + 2
        nsteps = get_printed_nsteps(traj_no)
        velocities = np.zeros((nsteps,n_atoms,3))
        c=0
        i=0
        for line in f:
            if c<2:
                c+=1
                continue
            elif c == lines_per_traj-1:
                i+=1
                c=0
                continue
            else:
                vel = np.array((float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])))
                velocities[i,c-2,:] = vel
                c+=1
                continue
    return velocities/(10*fs)

def parse_traj_positions(traj_no):
    """reads in fort.(1000+traj_no) and returns positions for all atoms and 
    time steps"""

    filename = 'fort.{}'.format(1000+traj_no)

    with open(filename) as f:
        n_atoms = get_natoms(traj_no)
        lines_per_traj = n_atoms + 2
        nsteps = get_printed_nsteps(traj_no)
        positions = np.zeros((nsteps,n_atoms,3))
        c=0
        i=0
        for line in f:
            if c<2:
                c+=1
                continue
            elif c == lines_per_traj-1:
                i+=1
                c=0
                continue
            else:
                pos = np.array((float(line.split()[1]),float(line.split()[2]),float(line.split()[3])))
                positions[i,c-2,:] = pos
                c+=1
                continue
    return positions

def get_element_list(traj_no):
    """reads in fort.(1000+traj_no) and returns list of atoms"""

    filename = 'fort.{}'.format(1000+traj_no)
    
    with open(filename) as f:
        n_atoms = get_natoms(traj_no)
        lines_per_traj = n_atoms + 2
        element_list = []
        c=0
        for line in f:
            if c<2:
                c+=1
                continue
            else:
                element_list.append(line.split()[0])
                if c == lines_per_traj-1:
                    break
                c+=1

    return element_list

def parse_unit_cell():
    filename = 'box_EFT'

    unit_cell = np.loadtxt(filename)



    return unit_cell

def get_printed_nsteps(traj_no):
    filename = 'fort.{}'.format(1000+traj_no)

    with open(filename) as f:
        n_atoms = get_natoms(traj_no)
        lines_per_traj = n_atoms + 2
        length_of_file = sum(1 for line in f)
        nsteps = int(length_of_file/lines_per_traj)


    return nsteps

def build_atoms_list(traj_no):

    positions = parse_traj_positions(traj_no)
    velocities = parse_traj_velocities(traj_no)
    element_list = get_element_list(traj_no)
    unit_cell = parse_unit_cell()
    nsteps = get_printed_nsteps(traj_no)

    atoms_list = []
    for t in range(nsteps):
        atoms = Atoms(symbols=element_list,positions=positions[t,:,:],cell=unit_cell,pbc=[True,True,True])
        atoms.set_velocities(velocities[t,:,:])
        atoms_list.append(atoms)
       
    return atoms_list

def parse_traj_tensors(traj_no):
    """reads in fort.(2000+traj_no) and returns positions for all atoms and 
    time steps"""
    
    filename = 'fort.{}'.format(2000+traj_no)
    nsteps = get_printed_nsteps(traj_no)
    
    with open(filename) as f:
        for i,line in enumerate(f):
            if i == 1:
                ndim = len(line.split())
                break
    with open(filename) as f:    
        lines_per_traj = ndim+1
        tensors = np.zeros((nsteps,ndim,ndim))
        c=0
        t=0
        for line in f:
            if c == 0:
                c+=1
                continue

            else:
                float_line = [float(i) for i in line.split()]
                tensors[t,c-1,:] = float_line
                if c == lines_per_traj-1:
                    c = 0
                    t+=1
                else:
                    c+=1
    # Tensor in g/mol/10fs 
    # g/mol = amu/atom
    return (tensors/(10*fs))

def mass_weight(tensor,friction_masses):

    ndim = np.shape(tensor)[0]

    for i in range(ndim):
        i_atom = i // 3
        for j in range(ndim):
            j_atom = j // 3
            tensor[i,j] = tensor[i,j]/np.sqrt(friction_masses[i_atom]*friction_masses[j_atom])

    return tensor

def project_tensor(tensor,atoms,friction_atoms):

    modes = calc_modes(atoms,friction_atoms)
    projected_tensor = np.dot(modes,np.dot(tensor,modes.transpose()))

    return projected_tensor

def parse_printed_time_step(traj_no):
    filename = 'fort.{}'.format(1000+traj_no)

    with open(filename) as f:
        n_atoms = get_natoms(traj_no)
        lines_per_traj = n_atoms + 2
        for i,line in enumerate(f):
            if i == 1:
                t1 = float(line.split()[1])
            elif i == 1+lines_per_traj:
                t2 = float(line.split()[1])
                break

        printed_time_step = t2 - t1

    return printed_time_step*10*fs

def parse_input_parameters():
    filename = glob.glob('*.inp')[0]
    print('Reading input parameters from:' + str(filename))

    with open(filename) as f:
        for line in f:
            #if 'Time' in line:

            if 'NNA' in line:
                Ni = int(line.split(',')[3])
                Ji = int(line.split(',')[4].split()[0])
                break

    return Ni,Ji

def bin_quantum(state):
    state = round(state)

    return state

def get_friction_masses(atoms,friction_atoms):

    friction_masses = atoms.get_masses()[friction_atoms]

    return friction_masses

def plot_relaxation_rates(fig,ax,traj_no,friction_atoms):
    """Plots relaxation rates of diagonal projected tensor elements
    across an entire trajectory. Saves to folder in cwd/N/J
    Also plots off_diagonals in separate file"""

    projected_tensors = get_mass_weighted_projected_tensors(traj_no,friction_atoms)    
    printed_time_step = parse_printed_time_step(traj_no)
    nsteps = get_printed_nsteps(traj_no)
    time_axis = np.arange(0,nsteps,1)*printed_time_step
    ndim = np.shape(projected_tensors)[1]

    for i in range(ndim):
        ax[0].plot(time_axis/ps,projected_tensors[:,i,i]*ps,label=labels[i],**line_settings)
        for j in range(i,ndim):
            if i == j:
                continue
            ax[1].plot(time_axis/ps,projected_tensors[:,i,j]*ps,label=labels[i]+labels[j],**line_settings)
    
    ax[0].set_ylim(0,2.5)
    ax[1].set_ylim(0,1.)
    plot_settings(ax[0])
    ax[0].tick_params(labelbottom=False)    
    plot_settings(ax[1])
    ax[1].legend(fontsize=10,fancybox=True,framealpha=0,loc=0,ncol=3)
    fig.text(-0.01, 0.5, r'Relaxation  rate / $\mathrm{ps}^{-1} $', va='center', rotation='vertical',fontsize=20)

    return

def fig_settings(fig):
    fig.set_figheight(22*0.393701)
    fig.set_figwidth(22*0.393701)
    return

def plot_settings(ax):
    # Hide the right and top spines
    fontsize=15

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.legend(fontsize=15,fancybox=True,framealpha=0,loc=0)
    ax.set_xlim(0)
    ax.set_ylim(0)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    return

def get_mass_weighted_projected_tensors(traj_no,friction_atoms):
    tensors = parse_traj_tensors(traj_no)
    atoms_list = build_atoms_list(traj_no)

    projected_tensors = np.zeros_like(tensors)
    friction_masses = get_friction_masses(atoms_list[0],friction_atoms)

    for i in range(np.shape(tensors)[0]):
        tensors[i,:,:]=mass_weight(tensors[i,:,:],friction_masses)
        projected_tensors[i,:,:] = project_tensor(tensors[i,:,:],atoms_list[i],friction_atoms)

    return projected_tensors

def get_projected_tensors(traj_no,friction_atoms):
    tensors = parse_traj_tensors(traj_no)
    atoms_list = build_atoms_list(traj_no)

    projected_tensors = np.zeros_like(tensors)

    for i in range(np.shape(tensors)[0]):
        projected_tensors[i,:,:] = project_tensor(tensors[i,:,:],atoms_list[i],friction_atoms)

    return projected_tensors

def calc_work(traj_no,friction_atoms):
    """Calculates the energy loss across trajectory
    for each mode"""
    ndim = len(friction_atoms)*3
    friction_projected_velocities = get_friction_projected_velocities(traj_no,friction_atoms)
    nsteps = get_printed_nsteps(traj_no)
    work = np.zeros((nsteps,ndim))


    friction_force_vecs = calculate_projected_forces(traj_no,friction_atoms)

    for i in range(nsteps):
        for j in range(ndim):
        #work[i] = np.dot(friction_projected_velocities[i,:],friction_force_vecs[i,:])
            work[i,j] = friction_projected_velocities[i,j]*friction_force_vecs[i,j]

    return work

def get_friction_projected_velocities(traj_no,friction_atoms):
    velocities = parse_traj_velocities(traj_no)
    ndim = len(friction_atoms)*3
    nsteps = get_printed_nsteps(traj_no)
    friction_projected_velocities = np.zeros((nsteps,ndim))
    atoms_list = build_atoms_list(traj_no)
    for i in range(nsteps):
        friction_projected_velocities[i,:] = velocity_transform(traj_no,atoms_list[i],velocities[i,friction_atoms,:],friction_atoms)

    return friction_projected_velocities

def velocity_transform(traj_no,atoms,velocities,friction_atoms):
    
    ndim = len(friction_atoms)*3
    velocity_tensor = np.zeros((ndim,ndim))

    for i in range(ndim):
        velocity_tensor[i,i]=velocities.flatten()[i]

    modes = calc_modes(atoms,friction_atoms)
    projected_tensor = np.dot(modes,np.dot(velocity_tensor,modes.transpose()))

    transformed_velocities = np.zeros((ndim))
    for i in range(ndim):
        transformed_velocities[i] = projected_tensor[i,i]
    return transformed_velocities

def calculate_projected_forces(traj_no,friction_atoms):
    ndim = len(friction_atoms)*3
    projected_tensors = get_projected_tensors(traj_no,friction_atoms)
    friction_projected_velocities = get_friction_projected_velocities(traj_no,friction_atoms)
    nsteps = get_printed_nsteps(traj_no)
    friction_force_vecs = np.zeros((nsteps,ndim))

    for i in range(nsteps):
        friction_force_vecs[i,:] = np.dot(projected_tensors[i,:,:],friction_projected_velocities[i,:])

    return friction_force_vecs

def plot_energy_loss(fig,ax,traj_no,friction_atoms):

    nsteps = get_printed_nsteps(traj_no)
    printed_time_step = parse_printed_time_step(traj_no)
    ndim = len(friction_atoms)*3
    work = calc_work(traj_no,friction_atoms)
    cumulative_work = np.zeros_like(work)
    total_work = np.zeros((nsteps))

    for j in range(ndim):
        cumulative_work[:,j] = np.cumsum(work[:,j])
        total_work[:] += cumulative_work[:,j]

    time_axis = np.arange(0,nsteps,1)*printed_time_step

    ax.plot(time_axis/ps,total_work,**{**line_settings,'marker':None},label='Total',color='black')
    for j in range(ndim):
        ax.plot(time_axis/ps,cumulative_work[:,j],**{**line_settings,'marker':None})#,label=labels[j])

    plot_settings(ax)
    ax.set_ylim(0,1.5)
    fig.text(0.47, 0.3, r'Energy loss / eV', va='center', rotation='vertical',fontsize=20)

    return

def plot_projected_velocities(fig,ax,traj_no,friction_atoms):

    friction_projected_velocities = get_friction_projected_velocities(traj_no,friction_atoms)

    nsteps = get_printed_nsteps(traj_no)
    printed_time_step = parse_printed_time_step(traj_no)
    ndim = len(friction_atoms)*3
    time_axis = np.arange(0,nsteps,1)*printed_time_step
 
    for j in range(ndim):
        ax.plot(time_axis/ps,friction_projected_velocities[:,j]*fs,label=labels[j],**{**line_settings,'marker':None})
    plot_settings(ax)
    ax.set_ylim(-0.04,0.04)
    fig.text(0.47, 0.7, r'Velocity / $\AA{}$ / fs', va='center', rotation='vertical',fontsize=20)

    return

def plot_traj_summary(traj_no,friction_atoms):
    print('Trajectory number = {}'.format(traj_no))
    instance_number = get_instance_number()
    fig, ax = plt.subplots(2, 2, sharex='all')#,constrained_layout=True)

    plot_relaxation_rates(fig,ax[:,0],traj_no,friction_atoms)
    plot_energy_loss(fig,ax[1,1],traj_no,friction_atoms)
    plot_projected_velocities(fig,ax[0,1],traj_no,friction_atoms)
    fig_settings(fig)
    ax[0,1].tick_params(labelbottom=False)
    ax[0,1].get_legend().remove()
    trapped = False
    try:
        lifetime,Nf,Jf,scat_angle = parse_traj_summary(traj_no)
    except:
        print('Trajectory number {} was not analysed in fort.26, it was trapped!'.format(traj_no))
        trapped = True
    
    Ni,Ji = parse_input_parameters()

    if not trapped:
        Nf = bin_quantum(Nf)
        Jf = bin_quantum(Jf)
        traj_text = r"""Lifetime = {:0.2f} fs, Scattering angle = {:0.2f}, N$_i$ = {}, N$_f$ = {}, J$_i$ = {}, J$_f$ = {}"""\
            .format(lifetime/fs,scat_angle,Ni,Nf,Ji,Jf)
        fig.text(0.5,0.92,traj_text,ha='center',fontsize=15)
        output_dir = 'Ni={}_Ji={}/Nf={}/Jf={}/'.format(str(Ni),str(Ji),str(Nf),str(Jf))
    else:
        fig.text(0.5,0.92,'Trapped',ha='center',fontsize=15)
        output_dir = 'Ni={}_Ji={}/trapped/'.format(str(Ni),str(Ji))
        ax[0,0].set_xlim(0,1)

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    fig.text(0.5, 0.03, "Time / ps", ha='center',fontsize=20)
    plt.subplots_adjust(wspace=0.5)
    fig.savefig(output_dir+'{:03d}_{:04d}_summary.pdf'.format(instance_number,traj_no),transparent=True,bbox_inches='tight')
    plt.close()
    return

#TODO: Write general information to file in order to calculate average stuff for each vib elastic, 1 quanta etc
#TODO: Only produce images for scattering angle <20 though they said didnt see a difference in the state to state scattering

def get_instance_number():

    filename = glob.glob('out*')[0]

    filename = filename.replace('out','')

    instance_number = int(filename)

    return instance_number
def calc_modes(atoms,friction_atoms):
    """Calculates required transformation matrix to convert diatomic
    friction tensor to internal coordinates as defined in Maurer et al PRL 2017
    """
    ndim = len(friction_atoms)*3
    modes = np.zeros([ndim,ndim])
    internals = np.zeros(6)
    f1 = friction_atoms[0]
    f2 = friction_atoms[1]

    pos1 = atoms.positions[f1]
    pos2 = atoms.positions[f2]

    mass1 = atoms.get_masses()[f1]
    mass2 = atoms.get_masses()[f2]
    
    com = (mass1*pos1[:] + mass2*pos2[:])/(mass1+mass2)



    vec = pos1[:] - pos2[:]
    norm = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]
    internals[0] = np.sqrt(norm)
    #phi
    internals[1] = np.arccos((vec[0]/np.sqrt(vec[0]*vec[0]+vec[1]*vec[1])))
    #theta
    vx_tilde = np.cos(internals[1])*vec[0]+np.sin(internals[1])*vec[1]
    internals[2] = np.arccos((vec[2]/np.sqrt(vec[2]*vec[2]+vx_tilde*vx_tilde)))
    internals[3:] = com[:]
    #mode 1 is the internal stretch
    modes[0,:3] = vec
    modes[0,3:] = -vec
    #mode 2 is 1st rotation angle
    # theta angle between molecular axis and z axis
    modes[1,:] = [-vec[1],vec[0],0.,vec[1],-vec[0],0.]
    #mode 3 is the 2nd rotation angle
    #it is defined as the crossing line between 2 planes
    #plane 1 is defined by the z-axis and the H2 axis
    #plane 2 is defined with the H2 axis as the normal vector
    vec2 = [0.,0.,0.]
    vec2[0] = np.cos(internals[1])*vec[2]
    vec2[1] = np.sin(internals[1])*vec[2]
    vec2[2] = -vx_tilde
    modes[2,:3] = np.array(vec2)
    modes[2,3:] = -np.array(vec2)

    #mode 4 is the x translation
    modes[3,:] = [1.,0.,0.,1.,0.,0.]
    #mode 5 is the y translation
    modes[4,:] = [0.,1.,0.,0.,1.,0.]
    #mode 6 is the z translation
    modes[5,:] = [0.,0.,1.,0.,0.,1.]

    for i in range(ndim):
        modes[i,:]/=np.linalg.norm(modes[i,:])
    return modes