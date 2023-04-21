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
from dscribe.descriptors import SOAP, ValleOganov, EwaldSumMatrix, LMBTR
from sklearn import preprocessing as skpreprocessing
from sklearn.utils.fixes import loguniform
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ExpSineSquared, WhiteKernel, DotProduct,  ConstantKernel, RationalQuadratic
from keras.models import Sequential
#from keras.layers import Dense, Activation, Dropout
from keras.models import Model   
from keras.layers import * 
from sklearn.model_selection import KFold
from tensorflow.keras.metrics import RootMeanSquaredError
from keras.optimizers import Adam, SGD, RMSprop
from tensorflow_addons.metrics import RSquare
from sklearn.decomposition import PCA
from joblib import dump
import glob

# import wandb
# from wandb.keras import WandbMetricsLogger, WandbCallback

import ase.db as adb
import pickle
import time
from ase.geometry import wrap_positions
from ase.io.trajectory import TrajectoryReader
plt.style.use('clb_publication')

descriptor_type = "invd"
use_scaler = True

def create_descriptor(atoms_list, mode="SOAP"):

    atoms_list = [atoms_list[i][[
    0,
    9,10,11,12,
    13,14,15,16,
    17,18,19,20,
    21,22,23,24
    ]] for i in range(len(atoms_list))]

    # atoms_list = [atoms_list[i] for i in range(len(atoms_list))]

    n_threads = 4

    elems = ['Pt', 'H']




    if mode == "SOAP":
        config={
            "nmax" : 3,
            "lmax" : 2,
            "rcut" : 6.,
            "average" : "off",
            "sigma" : 1.0
            }

        desc = SOAP(
            species = elems,    # List of the elements which the descriptor can accept
            periodic = True,    # We are only studying isolated molecules
            rcut = config["rcut"],         # Cutoff for local region in Angstroms (needs to be far enough to still include furthest hydrogens)
            nmax = config["nmax"],           # Number of radial basis functions
            lmax = config["lmax"],           # Maximum degree of spherical harmonics
            average=config["average"],       # If 'inner' averages over all sites to create a global descriptor.
            sigma = config["sigma"]
                        )

        X = create_soaps(atoms_list,desc, n_threads, config["average"], normalise=False, wrap_all=True)


    elif mode == "ValleOganov":
        vo = ValleOganov(
        species=elems,
        k2={
            "sigma": 10**(-0.5),
            "n": 100,
            "r_cut": 5
        },
        k3={
            "sigma": 10**(-0.5),
            "n": 100,
            "r_cut": 5
        },
        )

        X = vo.create(atoms_list)
    

    elif mode == "invd":

        X = create_inv_distance_matrix(atoms_list,wrap_all=True)


    else: 
        print("Unknown descriptor type")



    return X

def create_soaps(atoms_list,desc,n_threads,average='off', normalise=True,wrap_all=False):


    if wrap_all==True:
        [atoms_list[ii].wrap() for ii in range(len(atoms_list))]
    else:
        for ii in range(len(atoms_list)):
            atoms_list[ii].positions[0,:] = wrap_positions(atoms_list[ii].positions,atoms_list[ii].cell,pbc=True)[0,:]


    # h_idx = [[atoms.get_chemical_symbols().index('H')] for atoms in atoms_list]
    # soap_list = desc.create(atoms_list, positions=h_idx, n_jobs=n_threads)

    soap_list = desc.create(atoms_list, n_jobs=n_threads)

    # if i==0:
    if average=='off':
        X = soap_list[:,0,:]
    else:
        X = soap_list[:,:]
    # else:
        # if average=='off':
        #     X = np.append(X,soap_list[:,0,:],axis=0)
        # else:
        #     X = np.append(X,soap_list[:,:],axis=0)

    # Normalise SOAPs if requested.
    if normalise:
        scaler = skpreprocessing.StandardScaler().fit(X)
        X = scaler.transform(X)

    return X

def create_inv_distance_matrix(atoms_list,wrap_all=False):

    from ase.io.trajectory import TrajectoryReader

    dist = []
    inv_dist = []

    if wrap_all==True:
        [atoms_list[ii].wrap() for ii in range(len(atoms_list))]
    else:
        for ii in range(len(atoms_list)):
            atoms_list[ii].positions[0,:] = wrap_positions(atoms_list[ii].positions,atoms_list[ii].cell,pbc=True)[0,:]


    for j in range(len(atoms_list)):

        #dist is for inverse distance matrix later
        d = atoms_list[j].get_all_distances(mic=True)
        

        # d = d[[0,21,22,23,24],:]
        d = d[0,:]
        dist.append(d)
        
        natoms = len(atoms_list[j].get_masses())
        #np.fill_diagonal(d,1.0)
        d[0] = 1

        invd = np.ones((natoms,natoms)) / d
        np.fill_diagonal(invd,0.0)

        inv_dist.append(invd)


    inv_dist = np.array(inv_dist)
    inv_dist=inv_dist.reshape(len(dist),natoms*natoms)


    return inv_dist

def parse_frictions(filenames):

    for i,filename in enumerate(filenames):

        y1 = np.load(filename)

        if i == 0:
            y = y1
        else:
            y = np.append(y,y1,axis=0)

    return y

def kernel_ridge(X_train,y_train, X_test, y_test):

    model = KernelRidge(kernel="rbf",
        alpha=1e-1, gamma=1e-3)

    m = model.fit(X_train, y_train)

    y_pred = m.predict(X_test)
    y_pred_train = m.predict(X_train)

    root_mean_squared_error =  mean_squared_error(y_test, y_pred)
    r2_test =  model.score(X_test, y_test)

    print("RMSE: ", root_mean_squared_error)
    print("R2: ",r2_test)

    filename = 'KRRmodel.sav'
    pickle.dump(m, open(filename, 'wb'))


    return m

def grid_search_kernel(X_train,y_train,X_test,y_test):

    # model = GridSearchCV(
    #     MultiOutputRegressor(#
    #     KernelRidge()),
    #     param_grid={"estimator__kernel" : ['linear', 'poly','rbf'],
    #     "estimator__alpha": [1.,0.1, 1e-2, 1e-3,1e-4], 
    #     "estimator__gamma": np.logspace(-3, 3, 7),
    #     #"estimator__kernel__length_scale": np.logspace(2, 2,5),
    #     # "kernel__periodicity": np.logspace(0, 1,3),
    #     }, 
    #     cv=5)


    model = GridSearchCV(
        KernelRidge(),
        param_grid={"kernel" : ['linear', 'poly','rbf'],
        "alpha": [1.,0.1, 1e-2, 1e-3,1e-4], 
        "gamma": np.logspace(-3, 3, 7),
        #"estimator__kernel__length_scale": np.logspace(2, 2,5),
        # "kernel__periodicity": np.logspace(0, 1,3),
        }, 
        cv=5)

    m = model.fit(X_train, y_train)
    print(f"Best SVR with params: {m.best_params_} and R2 score: {m.best_score_:.3f}")


    y_pred = m.predict(X)
    np.save("friction_prune_predict.npy",y_pred)


    y_pred = m.predict(X_test)
    y_pred_train = m.predict(X_train)
    root_mean_squared_error =  mean_squared_error(y_test, y_pred)
    print("RMSE with test data: ps-1", root_mean_squared_error)

    filename = 'best_KRRmodel.sav'
    pickle.dump(m, open(filename, 'wb'))

    return m

def gpr(X_train,y_train,X_test,y_test):

    print("         ")
    print("--------------- GPR -------------------")

    # 1. 
    #kernel= 1. * Matern(length_scale_bounds=(1e4, 1e7))# + DotProduct() #+ RationalQuadratic() #+ WhiteKernel(noise_level_bounds=(1e-10,1e-1))
    #kernel= 1. * DotProduct() #+ WhiteKernel(noise_level_bounds=(1e-10,1e-1))

    # 2. bit smoother
    kernel= 1* RBF(length_scale_bounds=(1e-4, 1e3)) + WhiteKernel(noise_level_bounds=(1e-20,1e2)) 
    #kernel =  ConstantKernel(1e-20, (1e-25, 1e-15))* RBF(length_scale=1)
 
    #kernel=  1. *  DotProduct()#
    # 3. doesn't work, not positive semidefinite
    #kernel = 1.0 * ExpSineSquared(length_scale_bounds=(1e-3, 1e3), periodicity_bounds=(0., 1e3)) #+ WhiteKernel(1e-1)

    # 4. noisier (less smooth)
    #kernel= 1. * Matern(length_scale=1.0, length_scale_bounds=(1e-6, 1e1)) + DotProduct()

    # kernel = ConstantKernel(0.1, (1e-23, 1e10)) * \
    #      RBF(0.1*np.ones(X_train.shape[1]), (1e-10, 1e3) ) + \
    #      WhiteKernel(0.1, (1e-23, 1e5))

    model = GaussianProcessRegressor(kernel=kernel, alpha=1e-8, normalize_y=True,n_restarts_optimizer=2)
    start_time = time.time()
    m = model.fit(X_train, y_train)
    print(
        f"Time for GaussianProcessRegressor fitting: {time.time() - start_time:.3f} seconds"
    )

    y_pred = m.predict(X)
    np.save("friction_prune_predict.npy",y_pred)


    y_pred, y_pred_std = m.predict(X_test, return_std=True)
    np.save("friction_prune_std.npy", y_pred_std)

    r2_train = m.score(X_train,y_train)
    r2_test = m.score(X_test,y_test)


    print("         ")
    # print("R2 for training data ", str(r2_train))
    # print("R2 for test data", str(r2_test))

    #print(y_pred)

    y_pred_train = m.predict(X_train)

    print("RMSE with training data: ps-1", mean_squared_error(y_train, y_pred_train))
    print("RMSE with test data: ps-1", mean_squared_error(y_test, y_pred))


    print("Parameters: ",m.kernel_)
    filename = 'best_GPRmodel.sav'
    pickle.dump(m, open(filename, 'wb'))
    print("         ")
    print("----------------------------------")
    print("         ")
    return m

def mlpr(X_train, y_train, X_test, y_test):

    # Define the model
    model = MLPRegressor(max_iter=2000)

    # Define the hyperparameter grid
    param_grid = {
        'hidden_layer_sizes': [(10,), (50,), (100,)],
        'alpha': [0.001, 0.01, 0.1, 1.],
        'learning_rate_init': [0.001, 0.01, 0.1, 1.]
    }

    # Create the grid search object
    model = GridSearchCV(model, param_grid, cv=5, scoring='r2', n_jobs=-1)


    m = model.fit(X_train, y_train)
    print(f"Best SVR with params: {m.best_params_} and R2 score: {m.best_score_:.3f}")


    y_pred = m.predict(X)
    np.save("friction_prune_predict.npy",y_pred)


    y_pred = m.predict(X_test)
    y_pred_train = m.predict(X_train)
    root_mean_squared_error =  mean_squared_error(y_test, y_pred)
    print("RMSE with test data: ps-1", root_mean_squared_error)

    filename = 'best_mlprmodel.sav'
    pickle.dump(m, open(filename, 'wb'))

    return m

def parse_frictions(filenames):

    for i,filename in enumerate(filenames):

        y1 = np.load(filename)
        print("Parsed some frictions of length", str(len(y1)))

        if i == 0:
            y = y1
        else:
            y = np.append(y,y1,axis=0)

    return y

def plot_friction_traj(atoms,f,std=None,filename="test_traj_friction.pdf"):


    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    z = np.array([atoms[i].positions[0,2] for i in range(len(atoms))])

    dft = np.load("01_clean/test_traj/friction.npy")
    #dft = np.load("/Users/u1865573/work/2022/h_o_pt111/11_friction_ML/cluster/02_ox/test_traj/friction.npy")
    idx = [i for i in range(0,len(atoms),10)]

    skip = [110]
    #skip = [2,73]

    [idx.pop(i) for i in sorted(skip, reverse=True)]

    
    #indices = range(0,len(atoms),10)



    elements = np.shape(f)[1]

    # print(f[450:501,0])

    fig, ax = plt.subplots(elements, 1, figsize=(5, 6), sharex=True,sharey=True)

    ax2 = [ax[i].twinx() for i in range(elements)]

    idx2 = [i for i in range(0,len(atoms),1)]    

    for i in range(elements):

        ax[i].plot(idx2,f[:,i], linestyle='-', marker='None')
        ax2[i].plot(z-11.5, linestyle='--', color=colors[1], marker='None')

        ax[i].plot(idx,dft[:,i], linestyle="-", marker='None', mfc="none", markersize=1, color="black")

        if type(std) is np.ndarray:
            ax[i].fill_between(
            idx2,
            f[:,i]- 1.96* std[:,i],
            f[:,i] + 1.96 * std[:,i],
            #color="tab:green",
            color=colors[0],
            alpha=0.2,
            )




    #ax[-1].set_xlabel(r'$d$ / $\AA$')
    ax[-1].set_ylabel(r'Friction / ps$^{-1}$')
    ax2[-1].set_ylabel(r'$z$ / $\mathrm{\AA}$')


    #ax[0].set_xlim(0.,6.)
    ax[0].set_ylim(0,4.)

    fig.savefig(filename,transparent=False,bbox_inches='tight') 

def compare_to_dft(atoms,f,dft,tol=0.1, heights=[-10000.,10000.]):

    assert (np.shape(f)==np.shape(dft))
    print(np.shape(dft))

    bad = []
    for i in range(len(atoms)):

        z = atoms[i].positions[0,2]

        if z>heights[0] and z<heights[1]:

            for j in range(3):

                if abs(f[i,j] - dft[i,j])>tol:#(0.05*dft[i,j]):
                    bad.append(i)

    bad = np.unique(np.array(bad))

    bad_atoms = [atoms[i] for i in bad]

    return bad_atoms,f[bad,:], dft[bad,:]
                
def prune_f(f):
    " have to skip a few points and parse every 10th"

    idx = [i for i in range(0,len(f),10)]

    #skip = [110]
    skip = [2,73]

    [idx.pop(i) for i in sorted(skip, reverse=True)]

    return f[idx,:]

def prune_a(a):

    idx = [i for i in range(0,len(a),10)]

    #skip = [110]
    skip = [2,73]

    [idx.pop(i) for i in sorted(skip, reverse=True)]

    a_prune = [a[i] for i in idx]

    return a_prune

def write_atoms(atoms,calculation_status=2,filename="bad_atoms.db"):

    with adb.connect(filename) as f:

        for i in atoms:

            f.write(i,calculation_status=calculation_status)

def parse_atoms(filenames):

    from ase.io.trajectory import TrajectoryReader


    for i,db_file in enumerate(filenames):

        if ".db" in db_file:
            db = adb.connect(db_file)
            rows = db.select()
            atoms_list = [row.toatoms() for row in rows if(row.calculation_status==2)]

        else:
            traj = TrajectoryReader(db_file)

            atoms_list = [traj[t] for t in range(len(traj))]

            # atoms_list = []
            # for t in range(len(traj)):
            #     atoms_list.append(traj[t])


        if (i==0):
            atoms = atoms_list
        else:
            atoms+=atoms_list

        print("Parsed some Atoms of length", str(len(atoms_list)))


    return atoms

def plot_bad_atoms(atoms,f,dft,filename="bad_atoms.pdf"):

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    fig, ax = plt.subplots(3, 2, figsize=(5, 6), sharex=True,sharey=True, subplot_kw={'projection': '3d'})

    d1 = np.array([np.sort(atoms[i].get_all_distances(mic=True)[0,1:])[0] for i in range(len(atoms))])
    d2 = np.array([np.sort(atoms[i].get_all_distances(mic=True)[0,1:])[1] for i in range(len(atoms))])
    r =  np.array([atoms[i].positions[0,:] for i in range(len(atoms))])
    

    for i in range(3):

        
        # ax[i,0].plot(d1,f,linestyle="none", marker='o', color=colors[0])
        # ax[i,1].plot(d1,dft,color="black",linestyle="none", marker='x')

        # ax[i,0].scatter(d1,d2,f[:,i], marker='o', s=3)
        # ax[i,1].scatter(d1,d2,dft[:,i], marker='o', s=3) 

        ax[i,0].scatter(r[:,0], r[:,1], r[:,2]-11.5, c=f[:,i],marker='o', s=3, vmin=0, vmax=3)
        ax[i,1].scatter(r[:,0], r[:,1], r[:,2]-11.5, c=dft[:,i],  marker='o', s=3, vmin=0, vmax=3)


    # ax[0,0].set_xlim(0.8,3.5)
    # ax[0,0].set_ylim(0.8,3.5) 

    # plt.show()
    #fig.savefig(filename,transparent=False,bbox_inches='tight') 

def write_poorest_predictions(f, traj_filename, dft_filename, tol=1.0, heights=[-100000.,100000]):

    atoms = parse_atoms([traj_filename])
    # f = model.predict(X)

    #[atoms[i].wrap() for i in range(len(atoms))]
    for ii in range(len(atoms)):
        atoms[ii].positions[0,:] = wrap_positions(atoms[ii].positions,atoms[ii].cell,pbc=True)[0,:]

    dft = np.load(dft_filename)
    bad_atoms, bad_f, bad_dft = compare_to_dft(atoms,f,dft, tol=tol, heights=heights)
    print(str(len(bad_atoms)) + " bad predictions")

    write_atoms(bad_atoms,filename=traj_filename+"_bad.db")


    np.save(traj_filename+"_bad_friction_pred.npy", bad_f,)
    np.save(traj_filename+"_bad_friction_dft.npy", bad_dft)


    plot_bad_atoms(bad_atoms,bad_f,bad_dft, filename=traj_filename+"_bad_plot.pdf")

    return

def keras_model(X_train,y_train,X_test,y_test):

    config={
        "optimizer": 'rmsprop',
        "nodes":32,
        "epochs": 32,
        "learning_rate": 0.008031,
        "batch_size": 16,
        "regularization" : 0.002,
        "activation" : "relu",
        }

    learning_rate = config["learning_rate"]
    nodes = config["nodes"]
    optimizer = config["optimizer"]

    if config["optimizer"] == 'adam':
        optimizer = Adam(learning_rate=learning_rate)
    elif config["optimizer"] == 'sgd':
        optimizer = SGD(learning_rate=learning_rate)
    elif config["optimizer"] == 'rmsprop':
        optimizer = RMSprop(learning_rate=learning_rate)


    model = Sequential()

    model.add(Dense(units=nodes,input_shape=np.shape(X_train)[1:],activation=config["activation"])) 
    # model.add(BatchNormalization())
    model.add(Activation(config["activation"]))
    model.add(Dropout(config["regularization"]))

    model.add(Dense(units=nodes))
    # model.add(BatchNormalization())
    model.add(Activation(config["activation"]))
    model.add(Dropout(config["regularization"]))

    #model.add(Dense(nodes, activation=config["activation"]))
    model.add(Dense(units=1,activation=None))
    #model.add(Activation(config["activation"]))

    model.compile(
        #loss='binary_crossentropy',
        loss = 'mean_absolute_error', 
        optimizer=optimizer,             #also try adam
        metrics=[RSquare(),RootMeanSquaredError()]
        )
    history = model.fit(X_train, y_train, 
    epochs=config["epochs"],
    validation_split = 0.2,
    #validation_data=(X_test, y_test), 
    batch_size=config["batch_size"],) 
    # callbacks=[WandbMetricsLogger(),WandbCallback(model=model)])
    return model

def train(mode, X_train, y_train, X_test, y_test, grid_search=False):


    if mode=='KRR':

        if grid_search==True:
            model = grid_search_kernel(X_train,y_train,X_test,y_test)
        else:
            model = kernel_ridge(X_train,y_train,X_test,y_test)
        # KRR
        # filename = 'KRRmodel.sav'
        # model = pickle.load(open(filename, 'rb'))
        # y_pred = model.predict(X_traj)
        # plot_friction_traj(atoms,y_pred)

    elif mode=='GPR':
        # GPR
        model = gpr(X_train,y_train,X_test,y_test)

        filename = 'best_GPRmodel.sav'
        # model = pickle.load(open(filename, 'rb'))
        # print("Parameters: ", model.kernel_)

        # y_pred, y_pred_std = model.predict(X_traj, return_std=True)
        # plot_friction_traj(atoms,y_pred, std=y_pred_std)

    elif mode=='MLPR':
        # NN
        model = mlpr(X_train,y_train,X_test,y_test)
        filename = 'best_mlprmodel.sav'
        # model = pickle.load(open(filename, 'rb'))

        # y_pred = model.predict(X_traj)
        # plot_friction_traj(atoms,y_pred)

    elif mode=="keras":
        model = keras_model(X_train,y_train,X_test,y_test)
        # filename = 'keras.sav'
        # model = pickle.load(open(filename, 'rb'))
        #print("Mean/std X validation data", X.mean(), X.std())
        # y_pred = model.predict(X_traj)
        # plot_friction_traj(atoms,y_pred)



    return model

def plot_scatter_scaled(atoms_list,dft,f):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={"projection": "3d"})#, sharex=True,sharey=True)

    xf = np.array([atoms_list[i].get_scaled_positions(wrap=True)[0,0] for i in range(len(atoms_list))])
    yf = np.array([atoms_list[i].get_scaled_positions(wrap=True)[0,1] for i in range(len(atoms_list))]) 
    zf = np.array([atoms_list[i].get_scaled_positions(wrap=True)[0,2] for i in range(len(atoms_list))])


    xpt = np.array([atoms_list[0].get_scaled_positions(wrap=True)[-i,0] for i in range(1,5)])
    ypt = np.array([atoms_list[0].get_scaled_positions(wrap=True)[-i,1] for i in range(1,5)]) 
    zpt = np.array([atoms_list[0].get_scaled_positions(wrap=True)[-i,2] for i in range(1,5)])  


    # for i in range(np.shape(dft)[1]):
    #     ax[i].scatter(xf,yf,zf,c=dft[:,i])
        #ax[1,i].scatter(xf,yf,s=f[:,i])


    ax.scatter(xf,yf,zf,c=dft[:,0],vmin=1,vmax=2)
    ax.scatter(xpt,ypt,zpt,c="black",s=1000)
    plt.show()
    fig.savefig("scatter.pdf",transparent=False,bbox_inches='tight')


def plot_scatter(atoms_list,dft,filename="scatter.pdf"):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={"projection": "3d"})#, sharex=True,sharey=True)

    xf = np.array([atoms_list[i].positions[0,0] for i in range(len(atoms_list))])
    yf = np.array([atoms_list[i].positions[0,1] for i in range(len(atoms_list))]) 
    zf = np.array([atoms_list[i].positions[0,2] for i in range(len(atoms_list))])


    xpt = np.array([atoms_list[0].positions[-i,0] for i in range(1,5)])
    ypt = np.array([atoms_list[0].positions[-i,1] for i in range(1,5)]) 
    zpt = np.array([atoms_list[0].positions[-i,2] for i in range(1,5)])  


    # for i in range(np.shape(dft)[1]):
    #     ax[i].scatter(xf,yf,zf,c=dft[:,i])
        #ax[1,i].scatter(xf,yf,s=f[:,i])


    ax.scatter(xf,yf,zf,c=dft[:,0],vmin=1,vmax=2)
    ax.scatter(xpt,ypt,zpt,c="black",s=1000)
    plt.show()
    fig.savefig(filename,transparent=False,bbox_inches='tight')

def plot_more_friction_traj(atoms,f,std=None,filename="more_traj_friction.pdf"):


    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    z = np.array([atoms[i].positions[0,2] for i in range(len(atoms))])

    dft = np.load("01_clean/more_traj/friction.npy")



    root_mean_squared_error =  mean_squared_error(f, dft)
    print("RMSE more validation: ", root_mean_squared_error)

    for i in range(3):
        print("RMSE more validation" + str(i)+": ",  mean_squared_error(f[:,i], dft[:,i]))

    #dft = np.load("/Users/u1865573/work/2022/h_o_pt111/11_friction_ML/cluster/02_ox/test_traj/friction.npy")

    # skip = [110]
    # #skip = [2,73]

    # [idx.pop(i) for i in skip]

    
    #indices = range(0,len(atoms),10)



    elements = np.shape(f)[1]

    # print(f[450:501,0])

    fig, ax = plt.subplots(elements, 1, figsize=(14, 9), sharex=True,sharey=True)

    ax2 = [ax[i].twinx() for i in range(elements)]

    idx2 = [i for i in range(0,len(atoms),1)]    

    for i in range(elements):

        ax[i].plot(idx2,f[:,i], linestyle='-', marker='None')
        ax2[i].plot(idx2,z-11.5, linestyle='--', color=colors[1], marker='None')

        ax[i].plot(idx2,dft[:,i], linestyle="-", marker='None', mfc="none", markersize=1, color="black")

        if type(std) is np.ndarray:
            ax[i].fill_between(
            idx2,
            f[:,i]- 1.96* std[:,i], #* 1e5,
            f[:,i] + 1.96 * std[:,i],# * 1e5,
            #color="tab:green",
            color=colors[0],
            alpha=0.2,
            )




    #ax[-1].set_xlabel(r'$d$ / $\AA$')
    ax[-1].set_ylabel(r'Friction / ps$^{-1}$')
    ax2[-1].set_ylabel(r'$z$ / $\mathrm{\AA}$')


    #ax[0].set_xlim(0.,6.)
    ax[0].set_ylim(0,3.)

    fig.savefig(filename,transparent=False,bbox_inches='tight')

def plot_svd(X,y):



    pca=PCA(n_components=10)
    pca.fit(X)
    X_pca=pca.transform(X)


    fig, ax = plt.subplots(1, 1)
    ax.bar(np.arange(10),pca.explained_variance_ratio_)
    fig.savefig("pca.png", dpi=300, bbox_inches="tight", transparent=True)


    fig, ax = plt.subplots(3, 1, figsize=(8, 8))#, sharex=True,sharey=True)

    plt.xlabel("PC1")
    plt.ylabel("PC2")

    for i in range(3):
        ax[i].scatter(X_pca[:,0], X_pca[:,1], s=3, c=y[:,i], vmin=0., vmax=2.5, linewidths=0.4, label="datapoints")#, alpha=0.8)

    fig.savefig("pca_svd.pdf",transparent=False,bbox_inches='tight')

def plot_contour_prediction(model,scaler,atoms,z,filename="contour_prediction.pdf"):

    fig, ax = plt.subplots(4, 1, figsize=(4, 10))#, sharex=True,sharey=True)

    cell = atoms.cell

    frac = np.linspace(0,1,30)

    # z = 11.5

    print(cell)
    new = atoms

    b = atoms.positions

    traj = []
    pos = []
    for xf in frac:
        x = xf * cell[0,:]
        for yf in frac:
            y = yf * cell[1,:]

            xy = x + y

            pos.append(np.array([xy[0], xy[1], z]))


    


    y = np.zeros((len(pos),3))
    for i in range(len(pos)):
        new = atoms
        new.positions[0,:] = pos[i]

        X = create_descriptor([new,new], mode=descriptor_type)
        X = scaler.transform(X)

        for j in range(3):
            t = models[j].predict(X)
            y[i,j] = t[0]

    pos = np.array(pos)


    vmin = 0.
    vmax = 3.
    levels = np.linspace(vmin,vmax,51)
    for i in range(3):
        # y = models[i].predict(X)
        # print(y)
        tcf = ax[i].tricontourf(pos[:,0],pos[:,1],y[:,i],vmin=vmin,vmax=vmax,levels=levels,cmap="viridis",extend="both")
        tcf.cmap.set_over('white')
        tcf.cmap.set_under('navy')

    fig.colorbar(tcf, cax=ax[3], orientation='vertical')

    np.save("contour_pos.npy", pos)
    np.save("contour_pred.npy", y)
    fig.savefig(filename,transparent=False,bbox_inches='tight')
    
class multi_model():

    def __init__(self,models,cutoff):


        self.models = models
        self.cutoff = cutoff


    #How to employ cutoff in soap space

    def predict(self, descs):

        prediction = np.zeros((len(descs),len(self.models)))

        for i, desc in enumerate(descs):
            for m, model in enumerate(self.models):

                f = model.predict(desc.reshape(1,-1))

                f[f<0] = 0

                prediction[i,m] = f


        return prediction

def adaptive_sampling(filenames,models,scaler,outfile="highest_std_stuctures.db"):


    atoms = parse_atoms(filenames)

    X = create_descriptor(atoms, mode=descriptor_type)
    X = scaler.transform(X)



    f = np.zeros((len(atoms),3))
    std = np.zeros((len(atoms),3))

    for i in range(3):
        f[:,i], std[:,i] = models[i].predict(X, return_std=True)







    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    z = np.array([atoms[i].positions[0,2] for i in range(len(atoms))])

    elements = np.shape(f)[1]

    fig, ax = plt.subplots(elements, 1, figsize=(10, 8), sharex=True,sharey=True)

    ax2 = [ax[i].twinx() for i in range(elements)]

    idx2 = [i for i in range(0,len(atoms),1)]    

    for i in range(elements):

        ax[i].plot(idx2,f[:,i], linestyle='-', marker='None')
        ax2[i].plot(idx2,z-11.5, linestyle='--', color=colors[1], marker='None')

        if type(std) is np.ndarray:
            ax[i].fill_between(
            idx2,
            f[:,i]- 1.96* std[:,i] * 1e5,
            f[:,i] + 1.96 * std[:,i] * 1e5,
            #color="tab:green",
            color=colors[0],
            alpha=0.2,
            )

    #ax[-1].set_xlabel(r'$d$ / $\AA$')
    ax[-1].set_ylabel(r'Friction / ps$^{-1}$')
    ax2[-1].set_ylabel(r'$z$ / $\mathrm{\AA}$')


    #ax[0].set_xlim(0.,6.)
    ax[0].set_ylim(0,4.)

    fig.savefig('adaptive_trajs.pdf',transparent=False,bbox_inches='tight')


    bad_atoms = []

    for i in range(0,len(atoms),50):

        if np.any(std[i,:]*1e5 > 1.):

            bad_atoms.append(atoms[i])



    write_atoms(bad_atoms,calculation_status=0,filename=outfile)

def delete_boundary_structures(atoms,y):

    z = np.array([atoms[i].positions[0,2] for i in range(len(atoms))])
    d = np.array([np.sort(atoms[i].get_all_distances(mic=True)[0,:])[1] for i in range(len(atoms))]) 


    idx = np.argwhere((z>9.) & (z<21.5))[:,0]
    atoms = [atoms[i] for i in idx]
    y = y[idx]
    z = z[idx]

    for boundary in [0.0,0.5,1.0]:

        xf = np.array([atoms[i].get_scaled_positions(wrap=True)[0,0] for i in range(len(atoms))])
        idx_to_skip = np.argwhere(
            (np.isclose(xf,boundary,atol=0.01) == True) 
        )[:,0]

        y = np.delete(y,idx_to_skip, axis=0)
        [atoms.pop(i) for i in sorted(idx_to_skip, reverse=True)]


        yf = np.array([atoms[i].get_scaled_positions(wrap=True)[0,1] for i in range(len(atoms))])
        idx_to_skip = np.argwhere(
            (np.isclose(yf,boundary,atol=0.01) == True) 
        )[:,0]
        y = np.delete(y,idx_to_skip,axis=0)
        [atoms.pop(i) for i in sorted(idx_to_skip, reverse=True)]


    # SKip those where Pt atom goes higher
    idx_to_skip = np.argwhere(np.array([atoms[i].positions[1:,2] for i in range(len(atoms))])>11.8)[:,0]
    y = np.delete(y,idx_to_skip,axis=0)
    [atoms.pop(i) for i in sorted(idx_to_skip, reverse=True)]




        # idx4 = np.argwhere(np.isclose(yf[idx][idx2][idx3],0.5,atol=0.01) == True)

    return atoms, y

def plot_correlation(y_train,y_train_pred, y_test, y_test_pred):

    n_elements = np.shape(y)[1]

    fig, axs = plt.subplots(2, 1, figsize=(5, 10), sharex=True,sharey=True)

    for i in range(n_elements):

        axs[0].plot( y_train[:,i], y_train_pred[:,i], marker='x', linestyle='none', label=f'Element {i+1}')
        axs[1].plot(y_test[:,i], y_test_pred[:,i], marker='x', linestyle='none', label=f'Element {i+1}')

    axs[0].set_title('Training Data')
    axs[1].set_title('Test Data')


    perfect = np.linspace(0,10,1000)
    for ii,ax in enumerate(axs):
        if ii>0:
            ax.set_xlabel('Actual $\Lambda_{ii}$ / ps$^{-1}$')
        ax.set_ylabel('Predicted $\Lambda_{ii}$ / ps$^{-1}$')
        ax.legend()
        ax.set_xlim(0,np.ceil(np.max(y)))
        ax.set_ylim(0,np.ceil(np.max(y)))


        ax.plot(perfect,perfect,marker='None',color='black',linestyle='-',zorder=-6)

        
        # ax.annotate(soap_file.replace('_','-').replace('.npy',''),(2,0.1))
        #ax.annotate('hidden layers: '+str(hidden_layers[0]),(2.0,0.7))

        # for i in range(n_elements):
        #     ax.annotate('{} MAE: {:4.2f}'.format(i+1,maes[i, ii]),(3,3-i/2))
        #     ax.annotate('{} RMSE: {:4.2f}'.format(i+1,rmses[i, ii]),(3,1-i/2))

        fig.savefig('correlation.pdf',transparent=False,bbox_inches='tight')

def plot_pt_heights(atoms_list):

    heights = np.array([atoms[i].positions[9:,2] for i in range(len(atoms))])

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))#, sharex=True,sharey=True)


    

    for i in range(np.shape(heights)[1]):
        ax.hist(heights[:,i])
    

    fig.savefig("heights.pdf",transparent=False,bbox_inches='tight')





#friction_filenames = ['01_clean/friction.npy', "01_clean/more_traj/pes_database.db_bad_friction_dft.npy" ]
friction_filenames = [
"01_clean/friction.npy", 
"01_clean/more_training/friction.npy", 
"01_clean/fake_zero_structures/friction.npy",
"adaptive_sampling/01_clean/friction.npy"
]


y = parse_frictions(friction_filenames)

n_elements = np.shape(y)[1]


filenames = [
"01_clean/pes_database.db", 
"01_clean/more_training/pes_database.db", 
"01_clean/fake_zero_structures/zeros_database.db", 
"adaptive_sampling/01_clean/highest_std_stuctures.db"
]

atoms_list = parse_atoms(filenames)
[atoms_list[ii].wrap() for ii in range(len(atoms_list))]


print("Atoms list length", str(len(atoms_list)))
print("Friction length", str(len(y)))
assert(len(y)==len(atoms_list))


print('Total parsed data ' +str(len(atoms_list)))
atoms_list, y = delete_boundary_structures(atoms_list,y)
assert(len(y)==len(atoms_list))
print('Total training+test data ' +str(len(atoms_list)))

fig, ax = plt.subplots(3, 1, figsize=(7, 8))#, sharex=True,sharey=True)

z = np.array([atoms_list[i].positions[0,2] for i in range(len(atoms_list))])

for i in range(3):
    ax[i].plot(z, y[:,i], linestyle='none', markersize=0.1)
    ax[i].set_ylim(0,4.5)

fig.savefig("training.pdf",transparent=False,bbox_inches='tight')

X = create_descriptor(atoms_list, mode=descriptor_type)
print(np.shape(X))

# Split the data into training and validation sets
X_train, X_test, y_train, y_test, atoms_train, atoms_test = train_test_split(X, y, atoms_list, test_size=0.05)#, random_state=42)

print(len(X), len(X_train), len(atoms_train), len(X_test), len(atoms_test))

if use_scaler:
    # scaler = skpreprocessing.StandardScaler().fit(X_train)
    scaler =  skpreprocessing.MinMaxScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

print(np.max(X_train), np.max(X_test), np.min(X_train), np.min(X_test))
filename = 'xtrain_scaler.sav'
if use_scaler:
    pickle.dump(scaler, open(filename, 'wb'))

print("Mean/std X test data", X_test.mean(), X_test.std())

if use_scaler:
    X = scaler.transform(X)
plot_svd(X,y)

#     filename = 'best_KRRmodel_{}.sav'.format(["x","y","z"][i])
#     pickle.dump(m, open(filename, 'wb'))

traj_filenames = ["01_clean/test_traj/test_trajectory.traj"]


grid_search = False

# atoms_list = parse_atoms(traj_filenames)
atoms = TrajectoryReader(traj_filenames[0])
X_traj = create_descriptor(atoms, mode=descriptor_type)

if use_scaler:
    X_traj = scaler.transform(X_traj)

print("Mean/std X validation data", X_traj.mean(), X_traj.std())

y_pred = np.zeros((len(atoms),3))
y_pred_std = np.zeros((len(atoms),3))



mode = 'GPR'
atoms_more = parse_atoms(["01_clean/more_traj/pes_database.db"])
X_more = create_descriptor(atoms_more, mode=descriptor_type)
if use_scaler:
    X_more = scaler.transform(X_more)

y_more = np.zeros((len(atoms_more),3))
y_more_std = np.zeros((len(atoms_more),3))

y_train_pred = np.zeros((len(X_train),3))
y_test_pred = np.zeros((len(X_test),3))
y_both_pred = np.zeros((len(X),3))

models = []
for i in range(3):

    model = train(mode,X_train,y_train[:,i],X_test,y_test[:,i], grid_search=True)
    models.append(model)

    if mode=='GPR':
        y_pred_i, y_pred_std_i = model.predict(X_traj, return_std=True)
        y_pred_std[:,i] = y_pred_std_i

        y_train_pred_i = model.predict(X_train, return_std=False)
        y_test_pred_i = model.predict(X_test, return_std=False)
        y_both_pred_i = model.predict(X, return_std=False) 
    else:
        y_pred_i = model.predict(X_traj)
        y_train_pred_i = model.predict(X_train)
        y_test_pred_i = model.predict(X_test)
        y_both_pred_i = model.predict(X) 


    if mode == 'keras':
        y_pred[:,i] = y_pred_i[:,0]
        y_train_pred[:,i] =  y_train_pred_i[:,0]
        y_test_pred[:,i] =  y_test_pred_i[:,0]
        y_both_pred[:,i] = y_both_pred_i[:,0]
    else:
        y_pred[:,i] = y_pred_i
        y_train_pred[:,i] =  y_train_pred_i
        y_test_pred[:,i] =  y_test_pred_i
        y_both_pred[:,i] = y_both_pred_i


    # MORE

    if mode=='GPR':
        y_more_i, y_more_std_i = model.predict(X_more, return_std=True)
        y_more_std[:,i] = y_more_std_i
    else:
        y_more_i = model.predict(X_more)

    if mode == 'keras':
        y_more[:,i] = y_more_i[:,0]
    else:
        y_more[:,i] = y_more_i



plot_scatter(atoms_list,y)
plot_scatter(atoms_train,y_train,"scatter_pred_train.pdf")
plot_correlation(y_train, y_train_pred, y_test, y_test_pred)
plot_pt_heights(atoms_list)

fig, ax = plt.subplots(3, 1, figsize=(7, 8))#, sharex=True,sharey=True)
z = np.array([atoms_list[i].positions[0,2] for i in range(len(atoms_list))])
z_train = np.array([atoms_train[i].positions[0,2] for i in range(len(atoms_train))])
for i in range(3):
    ax[i].plot(z, y_both_pred[:,i], linestyle='none', mfc="none", markersize=0.1, color="red")
    ax[i].plot(z_train, y_train_pred[:,i], linestyle='none', mfc="none", markersize=0.1, color="black")
    ax[i].set_ylim(0,4.5)
fig.savefig("training_pred.pdf",transparent=False,bbox_inches='tight')



find_poorest = False
if find_poorest:
    write_poorest_predictions(y_more,
    traj_filename="01_clean/more_traj/pes_database.db",
    dft_filename="01_clean/more_traj/friction.npy",
    tol=0.2,
    heights=[9.,15.]
    )

    

if mode=='GPR':
    plot_friction_traj(atoms,y_pred, std=y_pred_std)
    plot_more_friction_traj(atoms_more,y_more,std=y_more_std)
    print("Standard deviation for ymore")
    print(y_more_std)
else:
    plot_friction_traj(atoms,y_pred)
    plot_more_friction_traj(atoms_more,y_more)





#  
# multimodel = multi_model(models,6.) 

# y_pred_multi = multimodel.predict(X_more)
# # plot_more_friction_traj(atoms_more,y_pred_multi)

# filename = mode+'_multimodel.sav'
# pickle.dump(multimodel, open(filename, 'wb'))


# dump(multimodel, mode+"_multimodel.lib")


for i,model in enumerate(models):
    filename = str(i)+mode+'_model.sav'
    pickle.dump(model, open(filename, 'wb')) 


# adapt_filenames = glob.glob("/Users/u1865573/work/2022/h_o_pt111/11_friction_ML/cluster/adaptive_sampling/scat*traj")
# adaptive_sampling(adapt_filenames,models,scaler)

plot_contour_prediction(models,scaler,atoms[0],10.,filename="contour_prediction_10.0.pdf")
plot_contour_prediction(models,scaler,atoms[0],11.5,filename="contour_prediction_11.5.pdf")
plot_contour_prediction(models,scaler,atoms[0],12.5,filename="contour_prediction_12.5.pdf")
plot_contour_prediction(models,scaler,atoms[0],13.5,filename="contour_prediction_13.5.pdf")



