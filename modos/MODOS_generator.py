import ase,os,shutil,sys
from ase.io import read,write
from ase.calculators.aims import Aims
from ase.visualize import view
import subprocess
from time import sleep

directory = sys.argv[1]
try:
    os.system("mkdir -p MODOS/"+directory)
except:
    print('MODOS already exists for this')
shutil.copy(directory+"/geometry.in", "MODOS/"+directory)
aims_srcMODOS = os.getenv("HOME")+'/fhi-aims.171221_1_release/src/MODOS/'

os.chdir(os.getenv("HOME")+"/MODOS_gen/MODOS/"+directory)
geom = read('geometry.in')
basis = 'tight'
kgrid = (4, 4, 1)
hpc = 'TINIS'
calc = Aims(aims_command = 'srun '+os.getenv("HOME")+'/software/fhiaims/code/binaries/'+hpc+'/aims.180128.scalapack.mpi.x',
            species_dir=os.getenv("HOME")+'/software/fhiaims/code/species_defaults/'+basis,
            xc='pbe',
            sc_accuracy_etot=1e-7,
            sc_accuracy_eev=1e-4,
            sc_accuracy_rho=1e-6,
            sc_accuracy_forces=1e-4,
            sc_iter_limit=100,
            spin='none',
            k_offset=(0.0, 0.0, 0.0),
            k_grid=kgrid,
            occupation_type = ['gaussian', 0.05],
            relativistic = 'atomic_zora scalar',
            mixer = 'pulay',
            KS_method = 'elpa',
            vdw_correction_hirshfeld = True,
            vdw_pair_ignore = ' Ag Ag',
            dos_kgrid_factors=kgrid,
            output=[' eigenvec_ovlp', ' dos -10.0 10.0 250 0.05'],
            use_full_spectrum= '.true.',
            collect_eigenvectors = '.true.',
            symmetry_reduced_k_grid = '.false.'
           )
geom.set_calculator(calc)
geom.get_potential_energy()
sleep(1)

try:
    os.mkdir(os.getenv("HOME")+'/MODOS_gen/MODOS/'+directory+'/overlayer/')
except:
    print('OVERLAYED FOLDER ALREADY MADE')
os.chdir(os.getenv("HOME")+'/MODOS_gen/MODOS/'+directory+'/overlayer/')
geom = read('../geometry.in')
del geom[[atom.index for atom in geom if atom.symbol=='Ag']]
geom.write('geometry.in')

geom = read('geometry.in')
basis = 'tight'
kgrid = (4, 4, 1)
hpc = 'TINIS'
calc = Aims(aims_command = 'srun '+os.getenv("HOME")+'/software/fhiaims/code/binaries/'+hpc+'/aims.180128.scalapack.mpi.x',
            species_dir=os.getenv("HOME")+'/software/fhiaims/code/species_defaults/'+basis,
            xc='pbe',
            sc_accuracy_etot=1e-7,
            sc_accuracy_eev=1e-4,
            sc_accuracy_rho=1e-6,
            sc_accuracy_forces=1e-4,
            sc_iter_limit=100,
            spin='none',
            k_offset=(0.0, 0.0, 0.0),
            k_grid=kgrid,
            occupation_type = ['gaussian', 0.05],
            relativistic = 'atomic_zora scalar',
            mixer = 'pulay',
            KS_method = 'elpa',
            dos_kgrid_factors=kgrid,
            output=[' eigenvec_ovlp', ' dos -10.0 10.0 250 0.05'],
            use_full_spectrum= '.true.',
            collect_eigenvectors = '.true.',
            symmetry_reduced_k_grid = '.false.'
           )
geom.set_calculator(calc)
geom.get_potential_energy()
sleep(1)

os.chdir(os.getenv("HOME")+"/MODOS_gen/MODOS/"+directory)
os.system(os.getenv("HOME")+'/fhi-aims.171221_1_release/src/MODOS/eigenvec_combine.out')
sleep(1)
print('Now doing the overlayer')
os.chdir(os.getenv("HOME")+"/MODOS_gen/MODOS/"+directory+'/overlayer')
os.system(os.getenv("HOME")+'/fhi-aims.171221_1_release/src/MODOS/eigenvec_combine.out')
os.system('mv eigenvec.out eigenvec.out_mol')
print('Overlayer complete')

try:
    os.mkdir(os.getenv("HOME")+"/MODOS_gen/MODOS/"+directory+'/datgen')
except:
    print('"datgen" folder already created previously. Moving there anyway')
os.chdir(os.getenv("HOME")+"/MODOS_gen/MODOS/"+directory+'/datgen')
os.system('cp ../eigenvec.out ../ovlp_mat.out ../overlayer/eigenvec.out_mol .')
os.system('cp '+aims_srcMODOS+'MODOS_control.in .')
os.system(aims_srcMODOS+'MODOS_FHI_aims_mpi')
print('All done, MODOS.dat files should be generated.')
