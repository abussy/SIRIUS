include("./Sirius.jl")
import .Sirius

using MKL
using BenchmarkTools
using MPI
using LinearAlgebra
@show BLAS.get_config()
BLAS.lbt_set_num_threads(2)

@show Sirius.libpath = ENV["LD_LIBRARY_PATH"]*"/libsirius.so"

MPI.Init()
@show comm = MPI.COMM_WORLD

@show Sirius.initialize(false)
@show ctx = Sirius.create_context_from_json(comm, "./sirius.json")
@show Sirius.initialize_context(ctx)


### Setting up the ground state SCF
@show k_grid, k_shift, use_symmetry = Sirius.get_kp_params_from_ctx(ctx)
@show kps = Sirius.create_kset_from_grid(ctx, k_grid, k_shift, use_symmetry)

@show gs = Sirius.create_ground_state(kps)

@show dtol, etol, itol, maxiter = Sirius.get_scf_params_from_ctx(ctx)
@show Sirius.find_ground_state(gs, true, true, dtol, etol, itol, maxiter)

@show temp, smear, kappa, tau, tol, maxiter, restart, pu = Sirius.get_nlcg_params_from_ctx(ctx)
@show converged = Sirius.nlcg(gs, kps, temp, smear, kappa, tau, tol, maxiter, restart, pu)

@show energy = Sirius.get_energy(gs, "total")

@show forces = Sirius.get_forces(gs, "total")

@show stress = Sirius.get_stress_tensor(gs, "total")

@show Sirius.free_ground_state_handler(gs)
@show Sirius.free_kpoint_set_handler(kps)
@show Sirius.free_context_handler(ctx)
@show Sirius.finalize(false)

MPI.Finalize()
