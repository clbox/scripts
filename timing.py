#Ideally record all timing outputs from aims output files and save into a easily plottable file
#Should also calculate average time for SCF, average time for CPSCF, average time for friction finite difference SCF and number
#of SCF steps and CPSCF steps
#Should work on multiple aims output files, recording also the number of atoms, tasks and k-points

import numpy as np
import sys




filenames = sys.argv[1:]
filenames.sort()

outfile = 'timing.out'

for filename in filenames:
    total_times = []
    headers = []

    n_scf = 0
    n_cpscf = 0
    n_tasks = 0
    n_atoms = 0
    n_scf_cycles = 0
    n_cpscf_cycles = 0
    use_dfpt = False
    #Whether to read line as timing output
    #Should only be active towards end of file
    read_total_times = False
    use_scalapack = False

    with open(filename,'r') as f:


        for line in f:

            if 'ELPA' in line:
                use_scalapack = True
                
            if 'Using' in line and 'parallel tasks' in line:
                n_tasks = int(line.split()[-3])
                continue

            if 'Number of atoms' in line:
                n_atoms = int(line.split()[-1])
                continue

            if '| Number of self-consistency cycles'  in line:
                n_scf = int(line.split()[-1])
                continue

            if 'Number of coupled perturbed self-consistency cycles' in line:
                n_cpscf = int(line.split()[-1])
                use_dfpt = True
                continue

            if '| Number of SCF (re)initializations' in line:
                n_scf_cycles = int(line.split()[-1])
                continue

            if '| Number of CPSCF (re)initializations' in line:
                n_cpscf_cycles = int(line.split()[-1])
                continue
            
            if 'Detailed time accounting' in line:
                read_total_times = True
                continue
            elif 'Partial memory accounting:' in line:
                read_total_times = False
                continue
            if read_total_times:
                if '|' in line:
                    total_times.append(float(line.split()[-4]))
                    headers.append(line.split('|')[1].split(':')[0].replace('  ','').replace(' ','_')[1:])

        with open(outfile,'a+') as f:
            f.write('Filename: '+filename+'\n')
            f.write('   Scalapack: '+str(use_scalapack)+'\n')
            f.write('   N_tasks: '+str(n_tasks)+'\n')
            f.write('   N_atoms: '+str(n_atoms)+'\n')
            f.write('   N_scf_iter: '+str(n_scf)+'\n')
            f.write('   N_scf_cycles: '+str(n_scf_cycles)+'\n')

            if use_dfpt:
                f.write('   DFPT in use \n')
                f.write('   N_cpscf_iter: '+str(n_cpscf)+'\n')
                f.write('   N_cpscf_cycles: '+str(n_cpscf_cycles)+'\n')

            for i,header in enumerate(headers):
                f.write('   '+header+'  '+str(total_times[i])+'\n')
