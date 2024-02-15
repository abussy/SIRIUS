include("SIRIUS.jl")
import .Sirius

using BenchmarkTools
using MPI
MPI.Init()
comm = MPI.COMM_WORLD

ctx = Sirius.context_handler(C_NULL)

@show Sirius.initialize(false)
@show Sirius.create_context_from_json(comm, ctx, "./sirius.json")
@show Sirius.initialize_context(ctx)


### Setting up the ground state SCF
kps = Sirius.kpoint_set_handler(C_NULL)
@show k_grid, k_shift, use_symmetry = Sirius.get_kp_info_from_ctx(ctx)
@show Sirius.create_kset_from_grid(ctx, k_grid, k_shift, use_symmetry, kps)

gs = Sirius.ground_state_handler(C_NULL)
@show Sirius.create_ground_state(kps, gs)

@show Sirius.find_ground_state(gs, true, true)#, density_tol=0.0001)

@show energy = Sirius.get_energy(gs, "total")

@show forces = Sirius.get_forces(gs, "total")

@show stress = Sirius.get_stress_tensor(gs, "total")

@show Sirius.free_ground_state_handler(gs)
@show Sirius.free_kpoint_set_handler(kps)
@show Sirius.free_context_handler(ctx)
@show Sirius.finalize(false)

MPI.Finalize()
