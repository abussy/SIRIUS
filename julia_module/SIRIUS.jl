#This is a tentative attempt at wrapping the SIRIUS C-style API into Julia
module Sirius

export initialize, finalize, create_context, context_handler, free_context_handler

using MPI
using JSON

#Hardcoded path to the SIRIUS library TODO: needs to be generalized
libpath = "/home/bussya/Documents/git/spack/opt/spack/linux-ubuntu23.10-skylake/gcc-13.2.0/sirius-develop-t5xeyb3ko6za5zyyjyd44ht4e36mgzxk/lib/libsirius.so"

#Strategy: we wrap each function of interest into a Julia function that handles optional and default arguments   

### Defining types to the handler pointers
mutable struct context_handler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct ground_state_handler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct kpoint_set_handler
   handler_ptr::Ref{Ptr{Cvoid}}
end

### Handler functions
function free_context_handler(ctx::context_handler)
   error_code__ = Ref{Cint}(0)
   #TODO: should we look into cfunctions? It looks like most people go that way
   @ccall libpath.sirius_free_object_handler(ctx.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_context_handler failed with error code", error_code__[])
   end
end

#TODO: we should probably add this in the struct as a destructor
function free_ground_state_handler(gs::ground_state_handler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(gs.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_ground_state_handler failed with error code", error_code__[])
   end
end

function free_kpoint_set_handler(kps::kpoint_set_handler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(kps.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_ground_state_handler failed with error code", error_code__[])
   end
end

function initialize_context(ctx::context_handler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_initialize_context(ctx.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.initialize_context failed with error code", error_code__[])
   end
end

### Simple init/finalize
function initialize(call_mpi_init::Bool)
   call_mpi_init__ = Ref{Cuchar}(call_mpi_init)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_initialize(call_mpi_init__::Ref{Cuchar}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.initialize failed with error code", error_code__[])
   end
end

function finalize(call_mpi_fin::Bool, call_device_reset::Bool=true, call_fftw_fin::Bool=true)
   call_mpi_fin__ = Ref{Cuchar}(call_mpi_fin)
   call_device_reset__ = Ref{Cuchar}(call_device_reset)
   call_fftw_fin__ = Ref{Cuchar}(call_fftw_fin)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_finalize(call_mpi_fin__::Ref{Cuchar}, call_device_reset__::Ref{Cuchar}, 
                                  call_fftw_fin__::Ref{Cuchar}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.finalize failed with error code", error_code__[])
   end
end

### Utility function that maps an integer to a C MPI_comm
#TODO: -is this complete? 
#      -should we allow C MPI_comms as argument in the C-style API to avoid back and forth?
comm2f(comm::MPI.Comm) = ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)

###TESTING a new function that initializes a context from a Json file
#TODO: make it without argument, and returning ctx?
function create_context_from_json(comm::MPI.Comm, ctx::context_handler, fname::String)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_context_from_json(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}},
                                                  fname::Cstring, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_context_from_json failed with error code", error_code__[])
   end
end

function get_kp_info_from_ctx(ctx::context_handler)

   k_grid = Vector{Cint}(undef, 3)
   k_shift = Vector{Cint}(undef, 3)
   use_symmetry = Ref{Cuchar}(false)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_kp_info_from_ctx(ctx.handler_ptr::Ptr{Cvoid}, k_grid::Ref{Cint}, 
                                              k_shift::Ref{Cint}, use_symmetry::Ref{Cuchar},
                                              error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_kp_info_from_ctx failed with error code", error_code__[])
   end

   return k_grid, k_shift, Bool(use_symmetry[])
end

#TODO: return a kpoins_set_handler rather than passing by arguments?
function create_kset_from_grid(ctx::context_handler, k_grid::Vector{Int32}, k_shift::Vector{Int32},
                               use_symmetry::Bool, kps::kpoint_set_handler)
   error_code__ = Ref{Cint}(0)
   use_symmetry__::Ref{Cuchar} = use_symmetry
   @ccall libpath.sirius_create_kset_from_grid(ctx.handler_ptr::Ptr{Cvoid}, k_grid::Ptr{Cint}, 
                                               k_shift::Ptr{Cint}, use_symmetry__::Ref{Cuchar},
                                               kps.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_kset_from_grid failed with error code", error_code__[])
   end
end

#TODO: return rather than pass by argument?
function create_ground_state(kps::kpoint_set_handler, gs::ground_state_handler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_ground_state(kps.handler_ptr::Ptr{Cvoid}, gs.handler_ptr::Ptr{Cvoid},
                                             error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_ground_state failed with error code", error_code__[])
   end
end

### Functions for the running of a SCF calculation
function find_ground_state(gs::ground_state_handler, initial_guess::Bool, save_state::Bool; #input
                           density_tol::Float64 = -1.0, energy_tol::Float64 = -1.0, #optional args 
                           iter_solver_tol::Float64 = -1.0, max_niter::Int64 = -1)  #if neg, take from ctx
   #input args
   initial_guess__ = Ref{Cuchar}(initial_guess)
   save_state__ = Ref{Cuchar}(save_state)
   
   #optional args
   if density_tol < 0.0
      density_tol__ = Ptr{Cdouble}(C_NULL)
   else
      density_tol__ = Ref{Cdouble}(density_tol)
   end 
   if energy_tol < 0.0
      energy_tol__ = Ptr{Cdouble}(C_NULL)
   else
      energy_tol__ = Ref{Cdouble}(energy_tol)
   end
   if iter_solver_tol < 0.0
      iter_solver_tol__ = Ptr{Cdouble}(C_NULL)
   else
      iter_solver_tol__ = Ref{Cdouble}(iter_solver_tol)
   end
   if max_niter < 0
      max_niter__ = Ptr{Cint}(C_NULL)
   else
      max_niter__ = Ref{Cint}(max_niter)
   end

   #output
   converged__ = Ref{Cuchar}(false)
   niter__ = Ref{Cint}(-1)
   rho_min__ = Ref{Cdouble}(-1.0)
   
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_find_ground_state(gs.handler_ptr::Ptr{Cvoid}, density_tol__::Ref{Cdouble},
                                           energy_tol__::Ref{Cdouble}, iter_solver_tol__::Ref{Cdouble},
                                           initial_guess__::Ref{Cuchar}, max_niter__::Ref{Cint},
                                           save_state__::Ref{Cuchar}, converged__::Ref{Cuchar}, 
                                           niter__::Ref{Cint}, rho_min__::Ref{Cdouble}, 
                                           error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.find_ground_state failed with error code", error_code__[])
   end
   
   return Bool(converged__[]), niter__[], rho_min__[]

end

function get_num_atoms(gs::ground_state_handler)
   num_atoms__ = Ref{Cint}(0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_num_atoms(gs.handler_ptr::Ptr{Cvoid}, num_atoms__::Ref{Cint},
                                       error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_num_atoms failed with error code", error_code__[])
   end
   return num_atoms__[]
end

function get_energy(gs::ground_state_handler, label::String)

   energy__ = Ref{Cdouble}(0.0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_energy(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, energy__::Ref{Cdouble},
                                    error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_energy failed with error code", error_code__[])
   end

   return energy__[]
end

function get_forces(gs::ground_state_handler, label::String)

   forces__ = Matrix{Cdouble}(undef, 3, get_num_atoms(gs))
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_forces(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, forces__::Ref{Cdouble},
                                    error_code__::Ref{Cint})::Cvoid 
   if error_code__[] != 0
      error("Sirius.get_forces failed with error code", error_code__[])
   end

   return forces__
end

function get_stress_tensor(gs::ground_state_handler,  label::String)
   
   stress__ = Matrix{Cdouble}(undef, 3, 3) 
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_stress_tensor(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, stress__::Ref{Cdouble},
                                           error_code__::Ref{Cint})::Cvoid 
   if error_code__[] != 0
      error("Sirius.get_stress_tensor failed with error code", error_code__[])
   end

   return stress__
end

### Below are a bunch of functions that were implemented for testing, but not really used in a SCF from JSON
#
#TODO: by default, comm_k and comm_band should be the same as comm? (need to verify behavior)
#      according to the Fortran API, if they are not explicitely passed, then they should be null pointers
function create_context(comm::MPI.Comm, ctx::context_handler, comm_k::MPI.Comm, comm_band::MPI.Comm)
   fcomm__::Int32 = comm2f(comm)
   fcomm_k__::Int32 = comm2f(comm_k)
   fcomm_band__::Int32 = comm2f(comm_band)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_context(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}}, 
                                        fcomm_k__::Int32, fcomm_band__::Int32, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_context failed with error code", error_code__[])
   end
end

function import_parameters(ctx::context_handler, parameters::String)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_import_parameters(ctx.handler_ptr::Ptr{Cvoid}, parameters::Cstring, 
                                           error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.import_parameters failed with error code", error_code__[])
   end
end

function set_lattice_vectors(ctx::context_handler, a1::Vector{Float64}, a2::Vector{Float64},
                             a3::Vector{Float64})
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_set_lattice_vectors(ctx.handler_ptr::Ptr{Cvoid}, a1::Ptr{Cdouble}, a2::Ptr{Cdouble}, 
                                             a3::Ptr{Cdouble}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.set_lattice_vectors failed with error code", error_code__[])
   end
end

function add_xc_functional(ctx::context_handler, name::String)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_add_xc_functional(ctx.handler_ptr::Ptr{Cvoid}, name::Cstring, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.add_xc_functional failed with error code", error_code__[])
   end
end


end #module
