# using Distributed, SlurmClusterManager
# addprocs(SlurmManager(); exeflags="--project")

# @everywhere begin
using NQCModels
using NQCBase  
using NQCDynamics
#using NNInterfaces
using NQCDynamics.InitialConditions: ConfigureAtomic
using FrictionProvidersLDFA
using NQCModels.AdiabaticModels
using NQCModels: FrictionModels
using PyCall
using EMTInterface
# end

using Unitful

using Plots
using Statistics: mean
using LinearAlgebra: norm, BLAS
using Libdl
using UnitfulAtomic
using Pandas: read_pickle

using SciMLBase
using Suppressor


BLAS.set_num_threads(1) # important for batching up onto nodes (embarssingly parallel)

# @pyimport ase.visualize as visualize

# using JLD2

# println("Workers ", nworkers())

# atoms = Atoms([:H, :H])
# cell = PeriodicCell([11.1175 -5.5588 0.0; 0.0 9.628 0.0; 0.0 0.0 70.3079])


function termination_condition(u, t, integrator)::Bool
    R = get_positions(u)
    zcom = au_to_ang(R[3,1])          # Convert vertical centre of mass to angstrom
    if zcom > 19.6                     # Scattering event
        return true
    elseif zcom < 0.
        return true
    else
        return false
    end
end

function calculate_ek(mass,vel)
    return 0.5.*sum(mass.*vel.*vel)
end

function write_trajectory(atoms,cell,traj; filename="test_trajectory.traj")

    traj_writer = io.trajectory.TrajectoryWriter(filename)

    #positions = zeros(size(R))

    for i_step=1:size(traj[:OutputPosition])[1]
        positions = traj[:OutputPosition][i_step]

        i_atoms = convert_to_ase_atoms(atoms,positions,cell)
        # io.write("i_atoms.in",i_atoms)

        traj_writer.write(i_atoms)
    end

end

function plot_theta_distribution(ensemble,filename="theta.png")

    # ....... ....     Plot scattering angle distribution
    nsamples = length(ensemble)
    θ_ss = zeros(n_scat)
    i = 1
    for s in 1:nsamples
        #pos = [ensemble[s].position[i][3,1] for i in 1:size(ensemble[s].position)[1]]
        # final_vel = ensemble[s].velocity[end][:,1]

        θ_s = ensemble[s][:OutputScatteringAngle]

        #final_z = ensemble[s][:OutputPosition][end][3,1]

        if ensemble[s][:OutputScatteredAtom]==1
            θ_ss[i] = θ_s
            i = i + 1
        end
        # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"Å",pos))
        # plot(auconvert.(u"fs",ensemble[s].t),auconvert.(u"eV",ensemble[s].kinetic))


    end
    histogram(θ_ss)
    xlabel!("θ_s")
    ylabel!("P(θ)")
    savefig(filename)


    θ_ss
end

function plot_final_kinetic_distribution(ensemble,atoms)

    # Plot height vs time for nsamples
    nsamples = length(ensemble)
    final_eks = zeros(nsamples)

    final_eks.=Inf
    # .... Histogram of kinetic energies
    mass_h=atoms.masses[1]
    for s in 1:nsamples
        final_vel = ensemble[s][:OutputVelocity][end][:,1]

        if ensemble[s][:OutputScatteredAtom]==true
            final_eks[s] = calculate_ek(mass_h,final_vel)
        end
    end


    b_range = range(-5, 5, length=71)
    histogram(auconvert.(u"eV",final_eks), bins=b_range)
    savefig("ek.png")


    final_eks
end

function plot_energy_loss_distribution(Ei, final_eks)

    b_range = range(-5, 5, length=71)
    histogram(ustrip.(Ei.-final_eks), bins=b_range)
    savefig("loss.png")
end

function calc_friction(fmodel,F,R)
    friction = NQCModels.friction!(fmodel,F,R)
    friction
end

function plot_friction_along_traj(fmodel,atoms,traj;filename="friction.png")

    f = zeros(size(traj[:Time])[1],3)
    z_pos = zeros(size(f))

    natoms = length(atoms.masses)
    mass_h=atoms.masses[1]



    for i_step=1:size(traj[:OutputPosition])[1]
        positions = traj[:OutputPosition][i_step]

        F = zeros(natoms*3,natoms*3)

        friction = calc_friction(fmodel,F,positions)

        for ii=1:3
            f[i_step,ii] = friction[ii,ii]
        end
        z_pos[i_step] = ustrip(auconvert(u"Å",positions[3,1]))
    end

    plot()

    for ii=1:3
        plot!(ustrip.(auconvert.(u"fs",traj[:Time])),
            ustrip.(auconvert.(u"ps^-1",f[:,ii])),
            legend=:topleft,label=ii)
    end

    # plot!(ensemble[i_step][:Time],z_pos)
    ylabel!("Λ / ps^-1")


    p = twinx()
    plot!(p,austrip.(traj[:Time]/u"fs"),z_pos,lc=:red,legend=:topright,label="R_z")
    xlabel!(p,"Time / fs")
    ylabel!(p,"Z / Å")

    savefig(filename)
end

function plot_friction_compare(fmodels,atoms,traj;filename="friction.png")

    f = zeros(size(traj[:Time])[1],3,length(fmodels))
    z_pos = zeros(size(f))

    natoms = length(atoms.masses)
    mass_h=atoms.masses[1]

    plot()


    for i_model in range(length(fmodels))
        for i_step=1:size(traj[:OutputPosition])[1]
            positions = traj[:OutputPosition][i_step]

            F = zeros(natoms*3,natoms*3)

            friction = calc_friction(fmodel,F,positions)

            for ii=1:3
                f[i_step,ii, i_model] = friction[ii,ii]
            end
            z_pos[i_step] = ustrip(auconvert(u"Å",positions[3,1]))
        end

        for ii=1:3
            plot!(ustrip.(auconvert.(u"fs",traj[:Time])),
                ustrip.(auconvert.(u"ps^-1",f[:,ii, i_model])),
                legend=:topleft,label=ii)
        end
    end

    # plot!(ensemble[i_step][:Time],z_pos)
    ylabel!("Λ / ps^-1")


    p = twinx()
    plot!(p,austrip.(traj[:Time]/u"fs"),z_pos,lc=:red,legend=:topright,label="R_z")
    xlabel!(p,"Time / fs")
    ylabel!(p,"Z / Å")

    savefig(filename)
end

function plot_friction_compare_height(filename="friction.pdf")

    emt_model, atoms, R, cell, frozen_atoms = set_up_clean_surface_LDFA()


    start = "/storage/chem/msrvhs/work/22/h_o_pt111/03_clean_surface_friction_convergence/01_hole_site/06/swap.in"
    odf_model_filename = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/11_friction_ML/best_GPRmodel.sav"
    model_scaler_filename = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/11_friction_ML/xtrain_scaler.sav"
    fmodel1 = set_up_odf_sklearn(start,odf_model_filename,model_scaler_filename)

    fmodel2 = set_up_pseudoldfa(emt_model,atoms)

    fmodels = [fmodel1 fmodel2]

    natoms = length(atoms.masses)


    heights = collect(LinRange(9,18,100))
    f = zeros(length(heights),3,length(fmodels))

    plot()
    for i_model in eachindex(fmodels)
        for i_step in eachindex(heights)

            R[3,1] = austrip(heights[i_step]*u"Å")
            F = zeros(natoms*3,natoms*3)
            friction = calc_friction(fmodels[i_model],F,R)

            for ii=1:3
                f[i_step,ii, i_model] = friction[ii,ii]
            end
        end

        for ii=1:3
            plot!(ustrip.((heights.-11.5)*u"Å"),
                ustrip.(auconvert.(u"ps^-1",f[:,ii, i_model])),
                legend=:topright,label=ii, lw=i_model)
        end
    end

    ylabel!("Λ / ps^-1")
    xlabel!("Height / Å")

    savefig(filename)
    library = Libdl.dlclose(emt_model.library) #Close Fortran code to hopefully reduce crashes due to lack
end

" Positions from Monte Carlo, velocities from Boltzmann distribution "
function initialize_from_montecarlo(v,r,chain,atoms,T,frozen, R, nsamples)

    natoms = length(atoms.masses)
    mobile = setdiff(2:natoms, frozen)

    # TODO: bottom layers set to zeros
    v_move = zeros(size(R[:,mobile]))
    v_froz = zeros(size(R[:,frozen]))

    vb = VelocityBoltzmann(T,atoms.masses[mobile],size(R[:,mobile]))

    temp_v = zeros(size(R))

    for i=1:nsamples
        v_move = rand(vb)
        rand_idx = rand(1:size(chain)[1])
        r[i] = hcat(r[i],chain[rand_idx][:,2:end])


        temp_v[:,1] =  v[i]
        temp_v[:,frozen] = v_froz
        temp_v[:,mobile] = v_move

        v[i] = temp_v

    end

    v,r
end

function set_up_clean_surface()

    slab = io.read("/storage/chem/msrvhs/work/22/h_o_pt111/03_clean_surface_friction_convergence/01_hole_site/06/swap.in")

    atoms, R, cell = convert_from_ase_atoms(slab)

    md_tian2_path = "/home/chem/msrvhs/git_repos/md_tian2/"
    lib_path=joinpath(md_tian2_path,"src/md_tian2_lib.so")
    pes_path = joinpath(md_tian2_path,"pes/EMT-HPt.pes")

    frozen_atoms = [collect(2:9);]

    model = NQCModels.md_tian2_EMT(atoms,cell, 
    lib_path,
    pes_path, 
    freeze=frozen_atoms
    )

    model,atoms,R,cell, frozen_atoms
end

function set_up_clean_surface_LDFA()

    #@pyimport ase.io as io

    io = pyimport("ase.io")

    slab = io.read("/storage/chem/msrvhs/work/22/h_o_pt111/03_clean_surface_friction_convergence/01_hole_site/06/swap.in")


    atoms, R, cell = convert_from_ase_atoms(slab)

    md_tian_path = "/home/chem/msrvhs/git_repos/md_tian2_LDFA/"
    lib_path=joinpath(md_tian_path,"md_tian_lib.so")

    md_tian2_path = "/home/chem/msrvhs/git_repos/md_tian2/" 
    pes_path = joinpath(md_tian2_path,"pes/EMT-HPt.pes")

    frozen_atoms = [collect(2:9);]

    # model = NQCModels.md_tian2_EMT(atoms,cell, 
    # lib_path,
    # pes_path, 
    # freeze=frozen_atoms
    # )

    model = EMTInterface.md_tian2_EMT(atoms,cell, 
    lib_path,
    pes_path, 
    freeze=frozen_atoms
    )

    return model,atoms,R,cell, frozen_atoms
end

function set_up_ox_surface()

    slab = io.read("/storage/chem/msrvhs/work/22/h_o_pt111/04_oxygen_surface_friction_convergence/01_hole_site/06/swap.in")

    atoms, R, cell = convert_from_ase_atoms(slab)
    
    md_tian_path = "/home/chem/msrvhs/git_repos/md_tian2_LDFA/"
    lib_path=joinpath(md_tian_path,"md_tian_lib.so")

    pes_path = joinpath(md_tian_path,"pes/PtOH.pes")
    
    frozen_atoms = [collect(3:10);]
    
    model = NQCModels.md_tian2_EMT(atoms,cell, 
    lib_path,
    pes_path, 
    freeze=frozen_atoms
    )
    
    model,atoms,R,cell, frozen_atoms
end

function set_up_odf_sklearn(initial_structure, model_filename, scaler_filename)
    dscr_d = pyimport("dscribe.descriptors")

    model_ml = read_pickle(model_filename)
    scaler = read_pickle(scaler_filename)
    ase_atoms = io.read(initial_structure)
    atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)

    println(ase_atoms)

    desc = dscr_d.SOAP(
        species = ["Pt", "O", "H"],
        periodic = true,
        rcut = 6,
        nmax = 3,
        lmax = 2,
        average="off" 
    )

    friction_model = SciKitFriction(desc, model_ml, ase_atoms, scaler; friction_unit=u"ps^-1")
    model = ODFriction(friction_model, atoms; friction_atoms=[1])

    model
end

function set_up_pseudoldfa(emt_model)

    fmodel = EMTFriction(emt_model;  friction_atoms=[1])

    return fmodel
end

function test_hcu_ldfa()
    slab = io.read("/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/12_hcu_benchmark/small_for_julia.in")

    atoms, R, cell = convert_from_ase_atoms(slab)

    cube_file = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/12_hcu_benchmark/cube_for_ldfa/cube_001_total_density.cube"
    density = CubeDensity(cube_file, cell, cell_matching_rtol=1e-3)

    frict_provider = LDFAFriction(density, atoms, friction_atoms=[1])

    natoms = length(atoms.masses)
    #F = zeros(natoms*3,natoms*3)

    heights = collect(LinRange(1,4,50))
    f = zeros(length(heights),3)



    plot()
    for i_step in eachindex(heights)

        R[3,1] = austrip((heights[i_step].+10.485)*u"Å")

        F = zeros(natoms*3,natoms*3)
        FrictionModels.friction!(frict_provider,F,R)

        for ii=1:3
            f[i_step,ii] = F[ii,ii]
        end
    end

    for ii=1:3
        plot!(ustrip.((heights)*u"Å"),
            ustrip.(auconvert.(u"ps^-1",f[:,ii])),
            legend=:topright,label=ii)
    end


    ylabel!("Λ / ps^-1")
    xlabel!("Height / Å")

    savefig("hcu_ldfa.pdf")


    #friction = NQCDynamics.Calculators.evaluate_friction!(sim.calculator,R)
end


function main()
    emt_model, atoms, R, cell, frozen_atoms = set_up_clean_surface_LDFA()
    # emt_model, atoms, R, cell, frozen_atoms = set_up_ox_surface()

    natoms = length(atoms.masses)


    # Friction model

        # 1. 
        # .............. ODF
        # start = "/storage/chem/msrvhs/work/22/h_o_pt111/03_clean_surface_friction_convergence/01_hole_site/06/swap.in"
        # odf_model_filename = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/11_friction_ML/best_GPRmodel.sav"
        # model_scaler_filename = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/11_friction_ML/xtrain_scaler.sav"
        # fmodel = set_up_odf_sklearn(start,odf_model_filename,model_scaler_filename)


        # 2.
        # ............ LDFA
       fmodel = set_up_pseudoldfa(emt_model)
        

    temperature = 293u"K"
    # sim = Simulation{Classical}(atoms, model; cell=cell, temperature=temperature)
    model = CompositeFrictionModel(emt_model,fmodel)
    sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=temperature)
    Ek = 1.92u"eV"
    z = 19.5u"Å"
    nsamples=100

    incidence_angle = 45.
    azimuthal_angle = 30.

    params = [incidence_angle,azimuthal_angle,Ek,temperature,nsamples]
    
    out_filename="HPt111_"
    for i in eachindex(params)
        out_filename *= string(ustrip(params[i]))*"_"
    end
    out_filename *= ARGS[1]
    out_filename *= ".h5"

    configurations = ConfigureAtomic.generate_configurations(sim, 
        samples=nsamples, translational_energy=Ek, height=z,incidence_angle=incidence_angle, azimuthal_angle=azimuthal_angle)

    v = first.(configurations)
    r = last.(configurations)


    # Monte Carlo
    # .... 1
    monte_sim = Simulation(atoms, emt_model; cell=cell, temperature=temperature)
    #Δ = Dict([(:H, 0.), (:Pt, 0.1)])
    Δ = Dict([(:H, 0.), (:O, 0.1), (:Pt, 0.1)])
    chain = InitialConditions.ThermalMonteCarlo.run_advancedmh_sampling(monte_sim, R, 1e3, Δ; move_ratio=0.9)

    # WRITE MONTE CARLo ------------------ 
    #traj_writer = io.trajectory.TrajectoryWriter("slab_trajectory.traj")

    # positions = zeros(size(R))
    # for i_step=1:size(chain)[1]
    # # for i_step=1:size(output.R)[1]
    #     positions = chain[i_step]
    #     i_atoms = convert_to_ase_atoms(atoms,positions,cell)
    #     traj_writer.write(i_atoms)
    # end

    v,r = initialize_from_montecarlo(v,r,chain,atoms,temperature,frozen_atoms, R, nsamples)

    #append starting surface to generated config of projectile atom
    #todo thermalize surface and freeze bottom atoms
    # v_surface = zeros(size(R))
    # for i=1:nsamples
    #     r[i] = hcat(r[i],R[:,2:end])
    #     v[i] = hcat(v[i],v_surface[:,2:end])
    # end


    distribution = DynamicalDistribution(v, r, size(R))


    timestep = 0.1u"fs"
    tspan = (0.0, 250.0u"fs")
    terminate = DynamicsUtils.TerminatingCallback(termination_condition)

    @suppress_err ensemble = run_dynamics(sim, tspan, distribution;selection=1:nsamples,
        dt=timestep, 
        output=(
            #OutputPosition, 
            #OutputVelocity,
            OutputPosition, 
            OutputVelocity,
            OutputScatteredAtom(19.6 * u"Å",1), 
            OutputScatteringAngle(sim; normal_vector=[0, 0, 1], incidence_angle=incidence_angle)),
        trajectories=nsamples, callback=terminate,
        reduction=FileReduction(out_filename))
        #ensemble_algorithm=SciMLBase.EnsembleDistributed())




    post_process = false


    if post_process
        println("----------- POST-PROCESSING ---------------")
        println(ensemble[1][:OutputScatteringAngle])


        # FIRST LOOP OVER ALL TRAj AND COUNT scattered

        n_scat = 0
        n_tot = 0
        for s in 1:nsamples
            n_scat = n_scat + ensemble[s][:OutputScatteredAtom]
            n_tot = n_tot + 1
        end
        n_trap = n_tot - n_scat

        println("Scattered: ", n_scat)
        println("Trapped: ", n_trap)



        final_eks = plot_final_kinetic_distribution(ensemble, atoms)

        loss = plot_energy_loss_distribution(Ek,auconvert.(u"eV",final_eks))

        i_traj = 1
        plot_friction_along_traj(fmodel,atoms,ensemble[i_traj]; filename="friction.png")

        # #  . ..... ..... Write Trajectory

        # for i_traj=1:nsamples
        i_traj = 1
            write_trajectory(atoms, cell, ensemble[i_traj]; filename="scat_traj_$i_traj.traj")
        # end

    end









    library = Libdl.dlclose(emt_model.library) #Close Fortran code to hopefully reduce crashes due to lack
    # # of deallocation statements
end

main()

# test_hcu_ldfa()
#plot_friction_compare_height()

