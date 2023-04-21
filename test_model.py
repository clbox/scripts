'''
Trains 3 Kernel Ridge Regression (KRR) models, one for each element
of the friction tensor.

Pretty much just uses the KernelRidge defaults.
More info available here: https://scikit-learn.org/stable/modules/generated/sklearn.kernel_ridge.KernelRidge.html#sklearn.kernel_ridge.KernelRidge

Might be a bit simpler with scikit-learn's MultiOutputRegressor,
which essentially does the same thing but neater: https://scikit-learn.org/stable/modules/multiclass.html#multioutput-regression
'''

import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
#from sklearn import preprocessing
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.kernel_ridge import KernelRidge
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.gaussian_process.kernels import Matern
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.multioutput import MultiOutputRegressor
from dscribe.descriptors import SOAP
from sklearn import preprocessing as skpreprocessing
import ase.db as adb
import pickle
from ase.io.trajectory import TrajectoryReader
plt.style.use('clb_publication')

def create_soaps(filenames,desc,n_threads,average='off', normalise=True):

    from ase.io.trajectory import TrajectoryReader

    for i,db_file in enumerate(filenames):

        if ".db" in db_file:
            db = adb.connect(db_file)
            rows = db.select()
            atoms_list = [row.toatoms() for row in rows if(row.calculation_status==2)]

        else:
            atoms_list = TrajectoryReader(db_file)


        h_idx = [[atoms.get_chemical_symbols().index('H')] for atoms in atoms_list]
        soap_list = desc.create(atoms_list, positions=h_idx, n_jobs=n_threads)

        if i==0:
            if average=='off':
                X = soap_list[:,0,:]
            else:
                X = soap_list[:,:]
        else:
            if average=='off':
                X = np.append(X,soap_list[:,0,:],axis=0)
            else:
                X = np.append(X,soap_list[:,:],axis=0)

        # Normalise SOAPs if requested.
        if normalise:
            scaler = skpreprocessing.StandardScaler().fit(X)
            X = scaler.transform(X)

    return X

def parse_frictions(filenames):

    for i,filename in enumerate(filenames):

        y1 = np.load(filename)

        if i == 0:
            y = y1
        else:
            y = np.append(y,y1,axis=0)

    return y


def plot_friction_traj(atoms,f,filename="test_traj_friction.pdf"):


    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    z = np.array([atoms[i].positions[0,2] for i in range(len(atoms))])

    elements = np.shape(f)[1]

    fig, ax = plt.subplots(elements, 1, figsize=(5, 6), sharex=True,sharey=True)

    ax2 = [ax[i].twinx() for i in range(elements)]

    for i in range(elements):

        ax[i].plot(f[:,i], linestyle='-', marker='none')
        ax2[i].plot(z-11.5, linestyle='--', color=colors[1], marker='none')




    #ax[-1].set_xlabel(r'$d$ / $\AA$')
    ax[-1].set_ylabel(r'Friction / ps$^{-1}$')
    ax2[-1].set_ylabel(r'$z$ / $\mathrm{\AA}$')


    #ax[0].set_xlim(0.,6.)
    ax[0].set_ylim(0,4.)

    fig.savefig(filename,transparent=False,bbox_inches='tight') 


config={
        "nmax" : 3,
        "lmax" : 2,
        "rcut" : 6.,
        "average" : "off",
        "sigma" : 1.0
        }
print(config)
elems = ['Pt', 'H']
desc = SOAP(
    species = elems,    # List of the elements which the descriptor can accept
    periodic = True,    # We are only studying isolated molecules
    rcut = config["rcut"],         # Cutoff for local region in Angstroms (needs to be far enough to still include furthest hydrogens)
    nmax = config["nmax"],           # Number of radial basis functions
    lmax = config["lmax"],           # Maximum degree of spherical harmonics
    average=config["average"],       # If 'inner' averages over all sites to create a global descriptor.
    sigma = config["sigma"])
                

n_threads = 4
#filenames = ['pes_database.db', '02_ox/pes_database.db']
#filenames = ['pes_database.db']
filenames = ["test_traj/test_trajectory.traj"]
X = create_soaps(filenames,desc, n_threads, config["average"])





#filename = 'best_KRRmodel_{}.sav'.format(["x","y","z"][i])

filename = 'best_KRRmodel.sav'

model = pickle.load(open(filename, 'rb'))

y_pred = model.predict(X)

atoms = TrajectoryReader(filenames[0])

plot_friction_traj(atoms,y_pred)