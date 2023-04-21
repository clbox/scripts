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
#using CubeLDFAModel
using FrictionProvidersLDFA

# atoms = Atoms([:H, :H])
# cell = PeriodicCell([11.1175 -5.5588 0.0; 0.0 9.628 0.0; 0.0 0.0 70.3079])
"Terminates simulation if returns `true`."
function termination_condition(u, t, integrator)::Bool
    R = get_positions(u)
    zcom = au_to_ang(R[3,1])          # Convert vertical centre of mass to angstrom
    if zcom > 8.00001                     # Scattering event
        return true
    else
        return false
    end
end

function calculate_ek(mass,vel)
    return 0.5.*sum(mass.*vel.*vel)
end

@pyimport ase.visualize as visualize
@pyimport ase.io as io

# slab = io.read("/storage/chem/msrvhs/work/22/julia_practice/h_pt111.in")
slab = io.read("/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/small.in")
#visualize.view(slab)
atoms, R, cell = convert_from_ase_atoms(slab)

cube_file = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/new_cube/cube_001_total_density.cube"
density = CubeDensity(cube_file, cell, cell_matching_rtol=1e-3)

# PES model
md_tian2_path = "/home/chem/msrvhs/git_repos/md_tian2/"
lib_path=joinpath(md_tian2_path,"src/md_tian2_lib.so")
pes_path = joinpath(md_tian2_path,"pes/EMT-HPt.pes")

pes_model = NQCModels.md_tian2_EMT(atoms,cell, 
lib_path,
pes_path
)


# .... LDFA model
#cube_file = "/storage/chem/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/nils_data/cube_001_total_density.cube"
# model = LDFAModel(model, cube_file, atoms, friction_atoms=[1], cell)

frict_provider = LDFAFriction(density, atoms, friction_atoms=[1])
model = CompositeFrictionModel(pes_model,frict_provider)

# .... test
# NQCModels.potential(model,R)
# NQCModels.derivative!(model,D,R)

sim = Simulation{Classical}(atoms, model; cell=cell)
Ek = 1.92u"eV"
z = 7.9u"Å"
nsamples=1000
configurations = ConfigureAtomic.generate_configurations(sim, 
    samples=nsamples, translational_energy=Ek, height=z)


v = first.(configurations)
r = last.(configurations)

#append starting surface to generated config of projectile atom
#todo thermalize surface and freeze bottom atoms
v_surface = zeros(size(R))
for i=1:nsamples
    r[i] = hcat(r[i],R[:,2:end])
    v[i] = hcat(v[i],v_surface[:,2:end])
end



distribution = DynamicalDistribution(v, r, (3,25)) #3 cart, 17 atoms
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")

timestep = 0.1u"fs"
tspan = (0.0, 250.0u"fs")
terminate = DynamicsUtils.TerminatingCallback(termination_condition)
# ensemble = run_ensemble(sim, tspan, distribution;selection=1:nsamples,
    # dt=timestep, output=(:position,:velocity,:kinetic), trajectories=nsamples, callback=terminate)
ensemble = run_dynamics(sim, tspan, distribution;selection=1:nsamples,
    dt=timestep, output=(OutputPosition, OutputVelocity), trajectories=nsamples, callback=terminate)

library = Libdl.dlclose(pes_model.library) #Close Fortran code to hopefully reduce crashes due to lack
# of deallocation statements

# Plot height vs time for nsamples
final_eks = zeros(nsamples)

# print(:OutputPosition)
# print(ensemble[1][:OutputPosition])



# .... Histogram of kinetic energies
mass_h=atoms.masses[1]
for s in 1:nsamples
    #pos = [ensemble[s].position[i][3,1] for i in 1:size(ensemble[s].position)[1]]
    # final_vel = ensemble[s].velocity[end][:,1]

    final_vel = ensemble[s][:OutputVelocity][end][:,1]
    final_eks[s] = calculate_ek(mass_h,final_vel)
    # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"Å",pos))
    # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"eV",ensemble[s].kinetic))


end
histogram(auconvert.(u"eV",final_eks))
savefig("test_ldfa.png")


# ...... Plot of z position
# for s in 1:nsamples
#     # x = ensemble[s][:Time]
#     # y = [ensemble[s][:OutputPosition][i][3,1] for i in 1:size(ensemble[s][:OutputPosition])[1]]
#     plot(ensemble[s][:Time],[ensemble[s][:OutputPosition][i][3,1] for i in 1:size(ensemble[s][:OutputPosition])[1]])
# end



# for s in 1:nsamples
#     # x = ensemble[s][:Time]
#     # y = [ensemble[s][:OutputPosition][i][3,1] for i in 1:size(ensemble[s][:OutputPosition])[1]]
#     plot!(ensemble[s][:Time],[ensemble[s][:OutputVelocity][i][3,1] for i in 1:size(ensemble[s][:OutputVelocity])[1]])
# end
# savefig("test_ldfa.png")

# plot(ensemble,:OutputPosition)

# ensemble[1].position[:]#[3,1]

