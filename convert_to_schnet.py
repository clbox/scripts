import numpy as np
from ase.db import connect
import schnetpack as spk
import sys
def string2array(string):
    """
    
    """
    dimension = len(string.split(']\n'))
    onedarray = np.fromstring((string.replace('[',' ').replace(']\n',' ')),dtype=float,sep=' ')
    return onedarray.reshape(dimension,dimension)



# PARSING DATA 
old_database = connect(sys.argv[1])

#available_properties = ["raw_coupling_energies","raw_coupling_elements","smear_frequency_energies","smeared_frequency_friction","markov_friction_tensor",]
available_properties = ["friction_tensor"]
dataset = spk.data.AtomsData("schnet.db", available_properties=available_properties)
raw_dataset = connect(old_database)

for idx in range(1,len(old_database)+1):
    print(idx)
    property_list = {}
    mol = raw_dataset.get_atoms(id=idx)
    #Markov
    row = raw_dataset.get(id=idx)
    try:
        markov_tensor = string2array(row.get('ft1'))
        property_list.update({"friction_tensor": markov_tensor})
    except:
        print('no markov')
        continue
    #RAW

    dataset.add_system(mol, **property_list)
    print('added')
