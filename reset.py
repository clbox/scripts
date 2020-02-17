from ase.io import read,write
from ase.visualize import view
import numpy as np
from ase import Atoms
import numpy as np
import os
from ase.db import connect
# coding: utf-8
from ase.calculators.aims import Aims
import os

con = connect('database.db')

for i in range(1,len(con)+1):
    row = con.get(id=i)   
    if row.cluster=='archer' or row.cluster=='archer2':
        con.update(i,calculation_status=0)

