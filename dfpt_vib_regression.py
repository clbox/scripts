import sys
import numpy as np

molecules = ['C2H2','C2H4','CH3Cl','CH4','H2CO','H2O2','N2H4','NH3','PH3','Si2H6','SiH4']
dimers = ['Cl2',  'ClF',  'CO',  'CS',  'F2',  'H2',  'HCl',  'HF',  'Li2',  'LiF',  'LiH',  'N2',  'Na2',  'NaCl',  'P2',  'SiO']
trimers = ['CO2',  'H2O',  'HCN',  'SH2',  'SO2']


reference_folder = sys.argv[1]
test_folder = sys.argv[2]

output_name = 'aims.out'
def parse_ir_data(directory,file):

    filename = directory+file
    read_data = False
    frequencies = []
    intensities = []
    with open(filename) as f:
        for line in f:
            if 'IR-intensity [D^2/Ang^2]' in line:
                read_data = True
                continue

            if read_data:

                if '------------------' in line or 'Finished' in line or len(line.split())<=0:
                    read_data = False
                    break
                frequencies.append(float(line.split()[1]))
                intensities.append(float(line.split()[2]))

    ir_data = np.zeros((len(frequencies),2))

    frequencies = np.array(frequencies)
    intensities = np.array(intensities)

    ir_data[:,0] = frequencies
    ir_data[:,1] = intensities

    return ir_data

print('1. Comparing molecules')

molecule_errors = np.zeros((len(molecules),2))
for i,molecule in enumerate(molecules):

    print(molecule)

    ref_ir_data = parse_ir_data(reference_folder,'molecules/'+molecule+'/output')

    test_ir_data = parse_ir_data(test_folder,'molecules/'+molecule+'/'+output_name)

    if np.shape(ref_ir_data)!=np.shape(test_ir_data):
        print('     *** ERROR - TEST DATA DIMENSIONS DOES NOT MATCH REF')
        print('     *** COULD BE FAILURE IN CALCULATON')
        continue

    frequency_mae = np.mean(np.abs(test_ir_data[:,0]-ref_ir_data[:,0]))
    intensity_mae = np.mean(np.abs(test_ir_data[:,1]-ref_ir_data[:,1]))

    print('     Frequency MAE / cm^(-1): ' + str(frequency_mae))
    print('     Intensity MAE / [D^2/Ang^2]: ' + str(intensity_mae))

    molecule_errors[i,0] = frequency_mae
    molecule_errors[i,1] = intensity_mae




print('2. Comparing dimers')

dimer_errors = np.zeros((len(dimers),2))
for i,dimer in enumerate(dimers):

    print(dimer)

    ref_ir_data = parse_ir_data(reference_folder,'dimer/'+dimer+'/output')

    test_ir_data = parse_ir_data(test_folder,'dimer/'+dimer+'/'+output_name)

    if np.shape(ref_ir_data)!=np.shape(test_ir_data):
        print('     *** ERROR - TEST DATA DIMENSIONS DOES NOT MATCH REF')
        print('     *** COULD BE FAILURE IN CALCULATON')
        continue


    frequency_mae = np.mean(np.abs(test_ir_data[:,0]-ref_ir_data[:,0]))
    intensity_mae = np.mean(np.abs(test_ir_data[:,1]-ref_ir_data[:,1]))

    print('     Frequency MAE / cm^(-1): ' + str(frequency_mae))
    print('     Intensity MAE / [D^2/Ang^2]: ' + str(intensity_mae))

    dimer_errors[i,0] = frequency_mae
    dimer_errors[i,1] = intensity_mae



print('3. Comparing trimers')

trimer_errors = np.zeros((len(trimers),2))
for i,trimer in enumerate(trimers):
    ref_ir_data = None
    test_ir_data = None

    print(trimer)

    ref_ir_data = parse_ir_data(reference_folder,'trimer/'+trimer+'/output')


    test_ir_data = parse_ir_data(test_folder,'trimer/'+trimer+'/'+output_name)

    if np.shape(ref_ir_data)!=np.shape(test_ir_data):
        print('     *** ERROR - TEST DATA DIMENSIONS DOES NOT MATCH REF')
        print('     *** COULD BE FAILURE IN CALCULATON')
        continue

    frequency_mae = np.mean(np.abs(test_ir_data[:,0]-ref_ir_data[:,0]))
    intensity_mae = np.mean(np.abs(test_ir_data[:,1]-ref_ir_data[:,1]))

    print('     Frequency MAE / cm^(-1): ' + str(frequency_mae))
    print('     Intensity MAE / [D^2/Ang^2]: ' + str(intensity_mae))

    trimer_errors[i,0] = frequency_mae
    trimer_errors[i,1] = intensity_mae


print('')
print('')
print('Frequency mean MAE for all molecules / cm^(-1):      '+str(np.mean(molecule_errors[:,0])))
print('Intensity mean MAE for all molecules / [D^2/Ang^2]:  '+str(np.mean(molecule_errors[:,1])))
print('')
print('Frequency mean MAE for all dimers / cm^(-1):         '+str(np.mean(dimer_errors[:,0])))
print('Intensity mean MAE for all dimers / [D^2/Ang^2]:     '+str(np.mean(dimer_errors[:,1])))
print('')
print('Frequency mean MAE for all trimers / cm^(-1):        '+str(np.mean(trimer_errors[:,0])))
print('Intensity mean MAE for all trimers / [D^2/Ang^2]:    '+str(np.mean(trimer_errors[:,1])))