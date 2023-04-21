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
using FrictionProvidersLDFA

# atoms = Atoms([:H, :H])
# cell = PeriodicCell([11.1175 -5.5588 0.0; 0.0 9.628 0.0; 0.0 0.0 70.3079])
"Terminates simulation if returns `true`."
function termination_condition(u, t, integrator)::Bool
    R = get_positions(u)
    zcom = au_to_ang(R[3,1])          # Convert vertical centre of mass to angstrom
    if zcom > 17.9                     # Scattering event
        return true
    elseif zcom < -1
        return true
    else
        return false
    end
end

function calculate_ek(mass,vel)
    return 0.5.*sum(mass.*vel.*vel)
end



function calc_friction(sim,R)
    friction = NQCDynamics.Calculators.evaluate_friction!(sim.calculator,R)
    friction
end


@pyimport ase.visualize as visualize
@pyimport ase.io as io

slab = io.read("/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/small.in")

atoms, R, cell = convert_from_ase_atoms(slab)

cube_file = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/new_cube/cube_001_total_density.cube"
density = CubeDensity(cube_file, cell, cell_matching_rtol=1e-3)

frict_provider = LDFAFriction(density, atoms, friction_atoms=[1])





md_tian2_path = "/home/chem/msrvhs/git_repos/md_tian2/"
lib_path=joinpath(md_tian2_path,"src/md_tian2_lib.so")
pes_path = joinpath(md_tian2_path,"pes/EMT-HPt.pes")

pes_model = NQCModels.md_tian2_EMT(atoms,cell, 
lib_path,
pes_path
)


model = CompositeFrictionModel(pes_model,frict_provider)
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=0u"K")
Ek = 1.92u"eV"
z = 17.8u"Å"
nsamples=10



configurations = ConfigureAtomic.generate_configurations(sim, 
    samples=nsamples, translational_energy=Ek, height=z,incidence_angle=45.)

v = first.(configurations)
r = last.(configurations)

#append starting surface to generated config of projectile atom
#todo thermalize surface and freeze bottom atoms
v_surface = zeros(size(R))
for i=1:nsamples
    r[i] = hcat(r[i],R[:,2:end])
    v[i] = hcat(v[i],v_surface[:,2:end])
end


distribution = DynamicalDistribution(v, r, size(R)) #3 cart, 17 atoms


timestep = 0.1u"fs"
tspan = (0.0, 250.0u"fs")
terminate = DynamicsUtils.TerminatingCallback(termination_condition)

ensemble = run_dynamics(sim, tspan, distribution;selection=1:nsamples,
    dt=timestep, 
    output=(
        OutputPosition, 
        OutputVelocity,
        OutputScatteredAtom(3.0u"Å",1), 
        OutputScatteringAngle(sim; normal_vector=[0, 0, 1])),
        # OutputDynamicsVariables),
    trajectories=nsamples, callback=terminate)

library = Libdl.dlclose(pes_model.library) #Close Fortran code to hopefully reduce crashes due to lack
# of deallocation statements



# FIRST LOOP OVER ALL TRAj AND COUNT scattered

n_scat = 0
n_tot = 0
for s in 1:nsamples
    global n_scat = n_scat + ensemble[s][:OutputScatteredAtom]
    global n_tot = n_tot + 1
end
n_trap = n_tot - n_scat

println("Scattered: ", n_scat)
println("Trapped: ", n_trap)


# Plot height vs time for nsamples
final_eks = zeros(n_scat)
# .... Histogram of kinetic energies
mass_h=atoms.masses[1]
i = 1
for s in 1:nsamples
    #pos = [ensemble[s].position[i][3,1] for i in 1:size(ensemble[s].position)[1]]
    # final_vel = ensemble[s].velocity[end][:,1]

    final_z = ensemble[s][:OutputPosition][end][3,1]
    final_vel = ensemble[s][:OutputVelocity][end][:,1]


    if ensemble[s][:OutputScatteredAtom]==1
        final_eks[i] = calculate_ek(mass_h,final_vel)
        global i = i + 1
    end
    # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"Å",pos))
    # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"eV",ensemble[s].kinetic))


end
histogram(auconvert.(u"eV",final_eks))
savefig("test.png")

#  . ..... ..... Write Trajectory


traj_writer = io.trajectory.TrajectoryWriter("test_trajectory.traj")

i_traj = 1
i_step = 1

positions = zeros(size(R))

for i_step=1:size(ensemble[i_traj][:OutputPosition])[1]
    positions = ensemble[i_traj][:OutputPosition][i_step]

    i_atoms = convert_to_ase_atoms(atoms,positions,cell)
    # io.write("i_atoms.in",i_atoms)

    traj_writer.write(i_atoms)
end

# ....... ....     Plot scattering angle distribution
θ_ss = zeros(n_scat)
i = 1
for s in 1:nsamples
    #pos = [ensemble[s].position[i][3,1] for i in 1:size(ensemble[s].position)[1]]
    # final_vel = ensemble[s].velocity[end][:,1]

    θ_s = ensemble[s][:OutputScatteringAngle]

    final_z = ensemble[s][:OutputPosition][end][3,1]

    if ensemble[s][:OutputScatteredAtom]==1
        θ_ss[i] = θ_s
        global i = i + 1
    end
    # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"Å",pos))
    # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"eV",ensemble[s].kinetic))


end
histogram(θ_ss)
xlabel!("θ_s")
ylabel!("P(θ)")
savefig("theta.png")




# ............ Recalculate friction and plot
i_traj = 1
i_step = 1

ldfas = zeros(size(ensemble[i_traj][:Time])[1])
z_pos = zeros(size(ldfas))

for i_step=1:size(ensemble[i_traj][:OutputPosition])[1]
    global positions = ensemble[i_traj][:OutputPosition][i_step]

    friction = calc_friction(sim,positions)
    ldfas[i_step] = friction[3,3]
    z_pos[i_step] = austrip(positions[3,1]/u"Å")
end

plot()
plot!(austrip.(ensemble[i_traj][:Time]/u"fs"),auconvert.(u"ps^-1",ldfas/mass_h)/u"ps^-1",legend=:topleft,label="LDFA")
# plot!(ensemble[i_step][:Time],z_pos)
ylabel!("LDFA / ps^-1")


p = twinx()
plot!(p,austrip.(ensemble[i_traj][:Time]/u"fs"),z_pos,lc=:red,legend=:topright,label="R_z")
xlabel!(p,"Time / fs")
ylabel!(p,"Z / Å")

savefig("ldfa.png")

