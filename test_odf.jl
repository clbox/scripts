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
@pyimport pickle
dscr_d = pyimport("dscribe.descriptors")

slab = io.read("/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/small.in")

atoms, R, cell = convert_from_ase_atoms(slab)

cube_file = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/07_nqcd_ldfa_benchmark/warwick_ldfa/new_cube/cube_001_total_density.cube"
density = CubeDensity(cube_file, cell, cell_matching_rtol=1e-3)

frict_provider = LDFAFriction(density, atoms, friction_atoms=[1])

model_filename = "/storage/mssgwp_grp/msrvhs/work/22/h_o_pt111/08_friction_ML/friction_ML/NNmodel.sav"




f = open(model_filename,"r")
odf_model = pickle.loads(pybytes(read(f)))


# GENERATE SOAP DESCRIPTER

desc = dscr_d.SOAP(
    species = ["Pt", "H"],
    periodic = true,
    rcut = 7.0,
    nmax = 8,
    lmax = 6,
    average="off" 
)




# friction = odf_model.predict(slab)


