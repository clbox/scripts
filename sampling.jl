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
using Distributions: Uniform


@pyimport ase.io as io

#slab = io.read("/storage/chem/msrvhs/work/22/julia_practice/h_pt111.in")
slab = io.read("/storage/chem/msrvhs/work/22/h_o_pt111/03_clean_surface_friction_convergence/01_hole_site/06/swap.in")
#visualize.view(slab)



atoms, R, cell = convert_from_ase_atoms(slab)

h_atom_index = 1

md_tian2_path = "/home/chem/msrvhs/git_repos/md_tian2/"
lib_path=joinpath(md_tian2_path,"src/md_tian2_lib.so")
pes_path = joinpath(md_tian2_path,"pes/EMT-HPt.pes")


frozen_atoms = [collect(2:9);]

model = NQCModels.md_tian2_EMT(atoms,cell, 
lib_path,
pes_path, 
freeze=frozen_atoms
)


temperature = 293u"K"


traj_writer = io.trajectory.TrajectoryWriter("slab_trajectory.traj")

# Monte Carlo
# .... 1
sim = Simulation(atoms, model; cell=cell, temperature=temperature)
Δ = Dict([(:H, 0.3), (:Pt, 0.1)])
R_new = R
positions = zeros(size(R))


# z_heights = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6] .* u"Å"
z_heights = range(8.5,17.5,15) .* u"Å"
z_heights = austrip.(z_heights)
for z in eachindex(z_heights)

    R_new[3,h_atom_index] = z_heights[z]
    lx = rand(Uniform(0.,1))
    ly = rand(Uniform(0.,1))

    vec1 = lx*cell.vectors[:,1]
    vec2 = ly*cell.vectors[:,2]
    vec = vec1 + vec2
    R_new[1,h_atom_index] = vec[1]
    R_new[2,h_atom_index] = vec[2]

    chain = InitialConditions.ThermalMonteCarlo.run_advancedmh_sampling(sim, R_new, 1e4, Δ; move_ratio=0.9)

    for i_step=1:size(chain)[1]
    # for i_step=1:size(output.R)[1]
        positions = chain[i_step]
        i_atoms = convert_to_ase_atoms(atoms,positions,cell)
        traj_writer.write(i_atoms)
    end
end




library = Libdl.dlclose(model.library) #Close Fortran code to hopefully reduce crashes due to lack
# of deallocation statements



