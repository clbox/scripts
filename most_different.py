#FROM ADAM MCSLOY
from scipy.spatial import distance_matrix
import numpy as np


def most_different(pool, count):
    """
    This uses the "Most Different Selection Method" to identify and return
    the most diverse selection of molecules possible from the available draw
    pool.

    Arguments:
        - pool (list): A list of ase.Atoms objects to select from.
        - count (int): The number of structures to return.

    Returns:
        - list: list of selected molecules
    """
    
    def get_most_different(draw_pool, selected_pool, dd_mat):
        # Identify the molecule most different to those already selected
        pl, sl, dif_mat = draw_pool, selected_pool, dd_mat

        # Get the index of the most different molecule and return it
        return pl.index(max(pl, key=lambda p: sum([dif_mat[p][s] for s in sl])))


    # Number of molecules to select must not be greater than the number supplied
    assert len(pool) >= count

    # Distance matrices of the molecules. The elements of these matrices are
    # the euclidean distances between the atoms.
    d_mats = [distance_matrix(i.positions, i.positions) for i in pool]

    # Difference matrix of different matrices. Saves us having to build these
    # repeatedly on the fly. (saves a lot of time)
    dd_mat = [[np.sum(np.square(i - j)) for j in d_mats] for i in d_mats]

    # Get the centroid of the distance matrices (i.e the average)
    cdm = np.average(d_mats, axis=0)

    # Identify molecule closest to this centroid (i.e. most average molecule)
    cent_index = np.argmin(np.array([np.sum(np.abs(d - cdm)) for d in d_mats]))

    # Create a pool from which to draw. Use integers to represent molecules so
    # that we can look them up in the pool, d_mats and dd_mats easily.
    draw_pool = list(range(len(pool)))

    # A list to hold the molecules that have been selected (integers of)
    selected_pool = []

    # Add centroid molecule to 'selected' list and remove it from the pool
    selected_pool.append(draw_pool.pop(cent_index))

    # Dummy loop, iteratively select the target number of molecules
    for _ in range(count - 1):

        # Find the molecule most different to those already selected
        diff_index = get_most_different(draw_pool, selected_pool, dd_mat)

        # Add to selected list and remove from pool list
        selected_pool.append(draw_pool.pop(diff_index))

    # Now pool all the selected molecules into a list and return it
    return [pool[i] for i in selected_pool]





if __name__ == '__main__':
    from ase.build import molecule


    # Build some dummy systems
    mols_in = []
    for _ in range(100):
        mol = molecule('H2O')
        mol.rattle(stdev=2)
        mols_in.append(mol)

    mols_out = most_different(mols_in, 10)

