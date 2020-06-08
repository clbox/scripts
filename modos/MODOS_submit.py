import sys,time,os
def subgen(name, nodes, code, hrs, mail):
    import os

    curdir = os.getcwd().split("/")
    if curdir[3] == "chem":
        dept = 'chem'
    if curdir[3] == "molases":
        dept = 'molases'
    if curdir[3] == "maths":
        dept = 'maths'
    if curdir[3] == 'diamond':
        dept = 'diamond'
    if curdir[1] == 'gpfs':
        procs = 28
        mem = 4571
        hpc = 'ORAC'
    if curdir[1] == 'tinisgpfs':
        procs = 16
        mem = 3882
        hpc = 'TINIS'

    module_list=['intel/2017.2.174-GCC-6.3.0-2.27 impi/2017.2.174  imkl/2017.2.174 \nulimit -s unlimited \nexport OMP_NUM_THREADS=1 \nexport MKL_NUM_THREADS=1 \nexport MKL_DYNAMIC=FALSE',
            'GCC/6.4.0-2.28 OpenMPI/2.1.1 \ntdir=$(mktemp -d /tmp/$SLURM_JOB_ID-XXXX) \n\ncp $SLURM_SUBMIT_DIR/*.inp $tdir/ \ncp $SLURM_SUBMIT_DIR/*.xyz $tdir/ \n\ncd $tdir',
            'GCC/7.3.0-2.30 OpenMPI/3.1.1 \ntdir=$(mktemp -d /tmp/$SLURM_JOB_ID-XXXX) \n\ncp $SLURM_SUBMIT_DIR/*.inp $tdir/ \ncp $SLURM_SUBMIT_DIR/*.xyz $tdir/ \n\ncd $tdir',
            'export OMP_NUM_THREADS=16 \n export PARNODES = 16' ]

    if code == 'aims':
          modules=module_list[0];
          path2bin="'srun /home/"+str(dept)+'/'+str(curdir[4])+'/software/fhiaims/code/binaries/'+str(hpc)+"/aims.180128.scalapack.mpi.x'"
          inp = ''
          post = ''    
    if code == 'MODOS':
          modules=module_list[0];
          path2bin="'python MODOS_generator.py "+name+"'"
          inp = ''
          post = ''
    if code == 'orca':
          modules=module_list[1]
          path2bin='/home/'+str(dept)+'/'+str(curdir[4])+'/software/ORCA/orca/orca'
          inp = 'input.inp'
          post = 'cp $tdir/* $SLURM_SUBMIT_DIR \n \ncd $SLURM_SUBMIT_DIR \nrm -rf $tdir'
    if code == 'orca412':
	  modules=module_list[1]
          path2bin='/home/'+str(dept)+'/'+str(curdir[4])+'/software/orca412/orca'
          inp = 'input.inp'
          post = 'cp $tdir/* $SLURM_SUBMIT_DIR \n \ncd $SLURM_SUBMIT_DIR \nrm -rf $tdir'
    if code == 'turbo':
          modules=module_list[3]
          path2bin = ''
          inp = ''
          post = ''

    sub = open('submit'+str(name)+'.sh', 'w')
    sub.write('#!/bin/bash \n'
        '#SBATCH -J '+str(name)+'\n'
        '#SBATCH --nodes='+str(nodes)+'\n'
        '#SBATCH --ntasks-per-node='+str(procs)+'\n'
        '#SBATCH --mem-per-cpu='+str(mem)+'\n'
        '#SBATCH --time='+str(hrs)+':00:00'+'\n'
        '#SBATCH --mail-type=END \n'
        '#SBATCH --mail-user='+str(mail)+'@warwick.ac.uk \n'
        '\n'
        'module load '+str(modules)+'\n'
        '\n'
        'export runbin='+str(path2bin)+'\n'
        '\n'
       '$runbin '+str(inp)+' > '+str(name)+'.out \n'
       '\n'+str(post)+'\n')
    sub.close()

subgen(sys.argv[1], sys.argv[2], 'MODOS', sys.argv[3], sys.argv[4])
print("All done, your file 'submit.sh' has been created!")
time.sleep(1)
os.system('sbatch submit'+str(sys.argv[1])+'.sh')
print('MODOS job running for '+str(sys.argv[1]))
