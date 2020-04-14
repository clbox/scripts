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

class venus_analysis():

    def __init__(self,traj_no,friction_atoms,mode=2):
        self.traj_no = traj_no
        self.friction_atoms = friction_atoms
        self.get_n_atoms()
        self.mode = mode

    def get_n_atoms(self):
        traj_no = self.traj_no
        filename = 'fort.{}'.format(1000+traj_no)

        with open(filename) as f:
            n_atoms = int(f.readline().split()[0])
        
        self.n_atoms = n_atoms

    def parse_traj_summary(self):
        """read fort.26 and gets all relevant trajectory information"""
        traj_no = self.traj_no
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

    def parse_traj_velocities(self,traj_no):
        """reads in fort.(1000+traj_no) and returns velocities for all atoms and 
        time steps"""
        filename = 'fort.{}'.format(1000+traj_no)

        with open(filename) as f:
            n_atoms = self.n_atoms
            lines_per_traj = n_atoms + 2
            nsteps = self.get_printed_nsteps()
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

    def parse_traj_positions(self,traj_no):
        """reads in fort.(1000+traj_no) and returns positions for all atoms and 
        time steps"""
        traj_no = self.traj_no
        filename = 'fort.{}'.format(1000+traj_no)

        with open(filename) as f:
            n_atoms = self.n_atoms
            lines_per_traj = n_atoms + 2
            nsteps = self.get_printed_nsteps()
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

    def get_element_list(self,traj_no):
        """reads in fort.(1000+traj_no) and returns list of atoms"""
        traj_no = self.traj_no
        filename = 'fort.{}'.format(1000+traj_no)
        
        with open(filename) as f:
            n_atoms = self.n_atoms
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

    def parse_unit_cell(self):
        filename = 'box_EFT'

        unit_cell = np.loadtxt(filename)



        return unit_cell

    def get_printed_nsteps(self):
        traj_no = self.traj_no
        filename = 'fort.{}'.format(1000+traj_no)

        with open(filename) as f:
            n_atoms = self.n_atoms
            lines_per_traj = n_atoms + 2
            length_of_file = sum(1 for line in f)
            nsteps = int(length_of_file/lines_per_traj)


        return nsteps

    def build_atoms_list(self):
        traj_no = self.traj_no
        positions = self.parse_traj_positions(traj_no)
        velocities = self.parse_traj_velocities(traj_no)
        element_list = self.get_element_list(traj_no)
        unit_cell = self.parse_unit_cell()
        nsteps = self.get_printed_nsteps()

        atoms_list = []
        for t in range(nsteps):
            atoms = Atoms(symbols=element_list,positions=positions[t,:,:],cell=unit_cell,pbc=[True,True,True])
            atoms.set_velocities(velocities[t,:,:])
            atoms_list.append(atoms)
        
        return atoms_list

    def parse_traj_tensors(self):
        """reads in fort.(2000+traj_no) and returns positions for all atoms and 
        time steps"""
        traj_no = self.traj_no
        filename = 'fort.{}'.format(2000+traj_no)
        nsteps = self.get_printed_nsteps()
        
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

    def mass_weight(self,tensor,friction_masses):
        
        ndim = np.shape(tensor)[0]

        for i in range(ndim):
            i_atom = i // 3
            for j in range(ndim):
                j_atom = j // 3
                tensor[i,j] = tensor[i,j]/np.sqrt(friction_masses[i_atom]*friction_masses[j_atom])

        return tensor

    def project_tensor(self,tensor,atoms):
        friction_atoms = self.friction_atoms
        modes = get_modes(atoms,friction_atoms,self.mode)
        projected_tensor = np.dot(modes,np.dot(tensor,modes.transpose()))

        return projected_tensor

    def parse_printed_time_step(self):
        traj_no = self.traj_no
        filename = 'fort.{}'.format(1000+traj_no)

        with open(filename) as f:
            n_atoms = self.n_atoms
            lines_per_traj = n_atoms + 2
            for i,line in enumerate(f):
                if i == 1:
                    t1 = float(line.split()[1])
                elif i == 1+lines_per_traj:
                    t2 = float(line.split()[1])
                    break

            printed_time_step = t2 - t1

        return printed_time_step*10*fs

    def parse_input_parameters(self):
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

    def bin_quantum(self,state):
        state = round(state)

        return state

    def get_friction_masses(self,atoms):
        friction_atoms = self.friction_atoms
        friction_masses = atoms.get_masses()[friction_atoms]

        return friction_masses

    def plot_relaxation_rates(self,fig,ax):
        """Plots relaxation rates of diagonal projected tensor elements
        across an entire trajectory. Saves to folder in cwd/N/J
        Also plots off_diagonals in separate file"""
        projected_tensors = self.get_mass_weighted_projected_tensors()    
        printed_time_step = self.parse_printed_time_step()
        nsteps = self.get_printed_nsteps()
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

    def get_mass_weighted_projected_tensors(self):
        tensors = self.parse_traj_tensors()
        atoms_list = self.build_atoms_list()

        projected_tensors = np.zeros_like(tensors)
        friction_masses = self.get_friction_masses(atoms_list[0])

        for i in range(np.shape(tensors)[0]):
            tensors[i,:,:]=self.mass_weight(tensors[i,:,:],friction_masses)
            projected_tensors[i,:,:] = self.project_tensor(tensors[i,:,:],atoms_list[i])

        return projected_tensors

    def get_projected_tensors(self):
        tensors = self.parse_traj_tensors()
        atoms_list = self.build_atoms_list()

        projected_tensors = np.zeros_like(tensors)

        for i in range(np.shape(tensors)[0]):
            projected_tensors[i,:,:] = self.project_tensor(tensors[i,:,:],atoms_list[i])

        return projected_tensors

    def calc_work(self):
        """Calculates the energy loss across trajectory
        for each mode"""
        friction_atoms = self.friction_atoms
        ndim = len(friction_atoms)*3
        friction_projected_velocities = self.get_friction_projected_velocities()
        nsteps = self.get_printed_nsteps()
        work = np.zeros((nsteps,ndim))


        friction_force_vecs = self.calculate_projected_forces()

        for i in range(nsteps):
            for j in range(ndim):
            #work[i] = np.dot(friction_projected_velocities[i,:],friction_force_vecs[i,:])
                work[i,j] = friction_projected_velocities[i,j]*friction_force_vecs[i,j]

        self.work = work
        return work

    def get_friction_projected_velocities(self):
        traj_no = self.traj_no
        friction_atoms = self.friction_atoms
        velocities = self.parse_traj_velocities(traj_no)
        ndim = len(friction_atoms)*3
        nsteps = self.get_printed_nsteps()
        friction_projected_velocities = np.zeros((nsteps,ndim))
        atoms_list = self.build_atoms_list()
        for i in range(nsteps):
            friction_projected_velocities[i,:] = self.velocity_transform(atoms_list[i],velocities[i,friction_atoms,:])

        return friction_projected_velocities

    def velocity_transform(self,atoms,velocities):
        friction_atoms = self.friction_atoms
        modes = get_modes(atoms,friction_atoms,self.mode)
        transformed_velocities = np.dot(modes,velocities.flatten())

        return transformed_velocities

    def calculate_projected_forces(self):
        friction_atoms = self.friction_atoms
        ndim = len(friction_atoms)*3
        projected_tensors = self.get_projected_tensors()
        friction_projected_velocities = self.get_friction_projected_velocities()
        nsteps = self.get_printed_nsteps()
        friction_force_vecs = np.zeros((nsteps,ndim))

        for i in range(nsteps):
            friction_force_vecs[i,:] = np.dot(projected_tensors[i,:,:],friction_projected_velocities[i,:])

        return friction_force_vecs

    def plot_energy_loss(self,fig,ax):
        friction_atoms = self.friction_atoms
        nsteps = self.get_printed_nsteps()
        printed_time_step = self.parse_printed_time_step()
        ndim = len(friction_atoms)*3
        work = self.calc_work()
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
        ax.set_ylim(0,2.2)
        fig.text(0.47, 0.3, r'Energy loss / eV', va='center', rotation='vertical',fontsize=20)

        return

    def plot_projected_velocities(self,fig,ax):

        friction_atoms = self.friction_atoms
        friction_projected_velocities = self.get_friction_projected_velocities()

        nsteps = self.get_printed_nsteps()
        printed_time_step = self.parse_printed_time_step()
        ndim = len(friction_atoms)*3
        time_axis = np.arange(0,nsteps,1)*printed_time_step
    
        for j in range(ndim):
            ax.plot(time_axis/ps,friction_projected_velocities[:,j]*fs,label=labels[j],**{**line_settings,'marker':None})
        plot_settings(ax)
        ax.set_ylim(-0.04,0.04)
        fig.text(0.47, 0.7, r'Velocity / $\AA{}$ fs$^{-1}$', va='center', rotation='vertical',fontsize=20)

        return

    def plot_traj_summary(self):
        traj_no = self.traj_no
        print('Trajectory number = {}'.format(traj_no))
        instance_number = self.get_instance_number()
        fig, ax = plt.subplots(2, 2, sharex='all')#,constrained_layout=True)

        self.plot_relaxation_rates(fig,ax[:,0])
        self.plot_energy_loss(fig,ax[1,1])
        self.plot_projected_velocities(fig,ax[0,1])
        fig_settings(fig)
        ax[0,1].tick_params(labelbottom=False)
        ax[0,1].get_legend().remove()
        self.trapped = False
        try:
            lifetime,Nf,Jf,scat_angle = self.parse_traj_summary()
        except:
            print('Trajectory number {} was not analysed in fort.26, it was trapped!'.format(traj_no))
            self.trapped = True
        
        Ni,Ji = self.parse_input_parameters()

        if not self.trapped:
            Nf = self.bin_quantum(Nf)
            Jf = self.bin_quantum(Jf)
            self.traj_text = r"""Lifetime = {:0.2f} fs, Scattering angle = {:0.2f}, N$_i$ = {}, N$_f$ = {}, J$_i$ = {}, J$_f$ = {}"""\
                .format(lifetime/fs,scat_angle,Ni,Nf,Ji,Jf)
            fig.text(0.5,0.92,self.traj_text,ha='center',fontsize=15)
            self.output_dir = 'Ni={}_Ji={}/Nf={}/Jf={}/'.format(str(Ni),str(Ji),str(Nf),str(Jf))
            self.summary_dir = 'Ni={}_Ji={}/Nf={}/'.format(str(Ni),str(Ji),str(Nf))
        else:
            fig.text(0.5,0.92,'Trapped',ha='center',fontsize=15)
            self.output_dir = 'Ni={}_Ji={}/trapped/'.format(str(Ni),str(Ji))
            self.summary_dir = self.output_dir
            ax[0,0].set_xlim(0,1)

        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        fig.text(0.5, 0.03, "Time / ps", ha='center',fontsize=20)
        plt.subplots_adjust(wspace=0.5)
        fig.savefig(self.output_dir+'{:03d}_{:04d}_summary.pdf'.format(instance_number,traj_no),transparent=True,bbox_inches='tight')
        plt.close()

        self.write_summary_to_file()

        return

    #TODO: Write general information to file in order to calculate average stuff for each vib elastic, 1 quanta etc
    #TODO: Only produce images for scattering angle <20 though they said didnt see a difference in the state to state scattering
    #TODO: Add incidence energy label to plot

    def get_instance_number(self):

        filename = glob.glob('out*')[0]

        filename = filename.replace('out','')

        instance_number = int(filename)

        self.instance_number = instance_number
        return instance_number

    def write_summary_to_file(self):
        traj_no = self.traj_no
        ndim = len(self.friction_atoms)*3
        total_work_dim = []
        total_work = 0
        fraction_energy_loss =[]
        for j in range(ndim):
            total_work_dim.append(np.sum(self.work[:,j]))
            total_work += np.sum(self.work[:,j])
        for j in range(ndim):
            fraction_energy_loss.append(total_work_dim[j]*100/total_work)
        with open(self.summary_dir+"summary.dat","a+") as f:
            f.write('Trajectory number = '+str(traj_no)+' Instance number = '+str(self.instance_number)+'\n')
            if not self.trapped:
                f.write(self.traj_text.replace('$','')+'\n')

            f.write('Total energy loss / eV, d = {:0.3f}, phi = {:0.3f}, theta = {:0.3f}, X = {:0.3f}, Y = {:0.3f}, Z = {:0.3f}, Total = {:0.3f}\n'\
                .format(*total_work_dim,total_work))

            f.write('Energy loss / %, d = {:0.1f}, phi = {:0.1f}, theta = {:0.1f}, X = {:0.1f}, Y = {:0.1f}, Z = {:0.1f}\n'\
                .format(*fraction_energy_loss))

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

def fig_settings(fig):
    fig.set_figheight(22*0.393701)
    fig.set_figwidth(22*0.393701)
    return

def get_modes(atoms,friction_atoms,mode=2):

    if mode == 1:
        modes = calc_modes(atoms,friction_atoms)

    elif mode == 2:
        modes = calc_modes2(atoms,friction_atoms)

    return modes
    
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

def calc_modes2(atoms,friction_atoms):
        """Calculates required transformation matrix to convert diatomic
        friction tensor to internal coordinates as defined by B. Jiang
        """

        ndim = len(friction_atoms)*3
        modes = np.zeros([ndim,ndim])
        f1 = friction_atoms[0]
        f2 = friction_atoms[1]

        pos1 = atoms.positions[f1]
        pos2 = atoms.positions[f2]

        x1 = pos1[0]
        y1 = pos1[1]
        z1 = pos1[2]
        x2 = pos2[0]
        y2 = pos2[1]
        z2 = pos2[2]

        m1 = atoms.get_masses()[f1]
        m2 = atoms.get_masses()[f2]

        mt = m1 + m2

        mr1 = m1/mt

        mr2 = m2/mt

        r = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        r1 = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        #mode 1 - r
        modes[0,0] = ((x1-x2)*mr2)/r
        modes[1,0] = ((y1-y2)*mr2)/r
        modes[2,0] = ((z1-z2)*mr2)/r

        modes[3,0] = ((x2-x1)*mr1)/r
        modes[4,0] = ((y2-y1)*mr1)/r
        modes[5,0] = ((z2-z1)*mr1)/r

        #mode 2  - theta
        modes[0,1] = ((x1-x2)*(z1-z2)*mr2)/r1
        modes[1,1] = ((y1-y2)*(z1-z2)*mr2)/r1
        modes[2,1] = -r1*mr2

        modes[3,1] = ((x2-x1)*(z1-z2)*mr1)/r1
        modes[4,1] = ((y2-y1)*(z1-z2)*mr1)/r1
        modes[5,1] = r1*mr1

        #mode 3 - phi
        modes[0,2] = -(y1-y2)*mr2
        modes[1,2] = (x1-x2)*mr2
        modes[2,2] = 0

        modes[3,2] = -(y2-y1)*mr1
        modes[4,2] = (x2-x1)*mr1
        modes[5,2] = 0

        #mode 4 is the x translation
        modes[:,3] = [1.,0.,0.,1.,0.,0.]
        #mode 5 is the y translation
        modes[:,4] = [0.,1.,0.,0.,1.,0.]
        #mode 6 is the z translation
        modes[:,5] = [0.,0.,1.,0.,0.,1.]

        return modes.transpose()