include("SIRIUS.jl")
import .Sirius

using MPI
MPI.Init()
comm = MPI.COMM_WORLD

ctx = Sirius.context_handler(C_NULL)

@show Sirius.initialize(false)
#@show Sirius.create_context(comm, ctx, comm, comm)
@show Sirius.create_context_from_json(comm, ctx, "./sirius.json")
@show Sirius.initialize_context(ctx)
#@show ctx
#@show Sirius.add_xc_functional(ctx, "PBE")
#@show Sirius.set_lattice_vectors(ctx, [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0])
#@show Sirius.import_parameters(ctx, "./sirius.json")

### Setting up the ground state SCF
kps = Sirius.kpoint_set_handler(C_NULL)
k_grid::Vector{Int32} = [2, 2, 2]
k_shift::Vector{Int32} = [1, 1, 1]
use_symmetry = false
@show Sirius.create_kset_from_grid(ctx, k_grid, k_shift, use_symmetry, kps)

gs = Sirius.ground_state_handler(C_NULL)
@show Sirius.create_ground_state(kps, gs)

@show Sirius.free_ground_state_handler(gs)
@show Sirius.free_kpoint_set_handler(kps)
@show Sirius.free_context_handler(ctx)
@show Sirius.finalize(false)

MPI.Finalize()
