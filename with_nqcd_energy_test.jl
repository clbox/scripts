# subroutine calc_force_wrapper(natoms, nbeads, ntypes, &
# r, f, epot, pes_file, natoms_list )
using NQCModels
using NQCBase
using PyCall
using NQCDynamics
#using NNInterfaces
using Unitful
using NQCDynamics.InitialConditions: ConfigureAtomic
using Plots
using Statistics: mean
using LinearAlgebra: norm
using Libdl
using UnitfulAtomic
# using JLD2
using CubeLDFAModel

#md_tian2 basic units:
#  Length : Ang
#  Time   : fs
#  Energy : eV




@pyimport ase.visualize as visualize
@pyimport numpy as np
@pyimport ase.io as io

slab = io.read("/storage/chem/msrvhs/work/22/h_o_pt111/06_nqcd_benchmark/mdtian2_test/bench/shift.poscar")



positions = slab.positions
atoms, R, cell = convert_from_ase_atoms(slab)

md_tian2_path = "/home/chem/msrvhs/git_repos/md_tian2/"
lib_path=joinpath(md_tian2_path,"src/md_tian2_lib.so")
pes_path = joinpath(md_tian2_path,"pes/EMT-HPt.pes")

model = NQCModels.md_tian2_EMT(atoms,cell, 
lib_path,
pes_path
)


print(potential(model,R))
print('\n')
austrip(0.14221007E+02 * u"eV")