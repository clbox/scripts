import numpy as np
from ase import Atoms
import schnetpack as spk

n_geoms = 50
available_properties = ["friction_tensor"]

dataset = spk.data.AtomsData("path_to_file.db", available_properties=available_properties)

for idx in range(n_geoms):
    geom = np.random.rand(5, 3)
    atypes = np.ones(5)

    mol = Atoms(atypes, geom)

    friction = np.random.rand(100, 21)

    properties = {
        "friction_tensor": friction
    }

    dataset.add_system(mol, **properties)
