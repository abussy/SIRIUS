#This is a tentative attempt at wrapping the SIRIUS C-style API into Julia
module Sirius

using MPI
using SIRIUS_jll

#lib = libsirius #from SIRIUS_jll
#lib = ENV["LD_LIBRARY_PATH"]*"/libsirius.so" #from spack local build
#@show lib
libpath = libsirius

### Utility function that maps an integer to a C MPI_comm
comm2f(comm::MPI.Comm) = ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)

### SIRIUS library inititialization and finalization
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

function is_initialized()
   status__ = Ref{Cuchar}(false)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_is_initialized(status__::Ref{Cuchar}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.is_initialized failed with error code", error_code__[])
   end
   return Bool(status__[])
end

### Defining types to the handler pointers of the C-API
mutable struct ContextHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct GroundStateHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct KpointSetHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

### Handler freeing function. Note: not added as finalizer as order metters
#TODO: could pass these destructors as struct finalizer, except that they may be triggered in any order
function free_context_handler(ctx::ContextHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(ctx.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_context_handler failed with error code", error_code__[])
   end
end

function free_ground_state_handler(gs::GroundStateHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(gs.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_ground_state_handler failed with error code", error_code__[])
   end
end

function free_kpoint_set_handler(kps::KpointSetHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(kps.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_kpoint_set_handler failed with error code", error_code__[])
   end
end

### Simulation context related functions
function initialize_context(ctx::ContextHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_initialize_context(ctx.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.initialize_context failed with error code", error_code__[])
   end
end

function create_context_from_json(comm::MPI.Comm, fname::String)
   ctx = ContextHandler(C_NULL)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_context_from_json(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}},
                                              fname::Cstring, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end

function create_context_from_json(fname::String)
   ctx = ContextHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_context_from_json_commworld(ctx.handler_ptr::Ref{Ptr{Cvoid}},
                                                        fname::Cstring, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end

function create_context(comm:: MPI.Comm)
   ctx = ContextHandler(C_NULL)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_empty_context(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}}, 
                                              error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_empty_context failed with error code", error_code__[])
   end
   return ctx
end

function get_kp_info_from_ctx(ctx::ContextHandler)

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

### Kpoint set related functions
function create_kset_from_grid(ctx::ContextHandler, k_grid::Vector{Int32}, k_shift::Vector{Int32},
                               use_symmetry::Bool)
   kps = KpointSetHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   use_symmetry__::Ref{Cuchar} = use_symmetry
   @ccall libpath.sirius_create_kset_from_grid(ctx.handler_ptr::Ptr{Cvoid}, k_grid::Ptr{Cint}, 
                                           k_shift::Ptr{Cint}, use_symmetry__::Ref{Cuchar},
                                           kps.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_kset_from_grid failed with error code", error_code__[])
   end
   return kps
end

### Ground state related functions
function create_ground_state(kps::KpointSetHandler)
   gs = GroundStateHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_ground_state(kps.handler_ptr::Ptr{Cvoid}, gs.handler_ptr::Ptr{Cvoid},
                                         error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_ground_state failed with error code", error_code__[])
   end
   return gs
end

function find_ground_state(gs::GroundStateHandler, initial_guess::Bool, save_state::Bool;
                           density_tol::Float64 = -1.0, energy_tol::Float64 = -1.0, 
                           iter_solver_tol::Float64 = -1.0, max_niter::Int64 = -1)

   # Optional keyword arguments: default behavior is to take value from the JSON file

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

###  Function for querying the results
function get_num_atoms(gs::GroundStateHandler)
   num_atoms__ = Ref{Cint}(0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_num_atoms(gs.handler_ptr::Ptr{Cvoid}, num_atoms__::Ref{Cint},
                                   error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_num_atoms failed with error code", error_code__[])
   end
   return num_atoms__[]
end

function get_energy(gs::GroundStateHandler, label::String)

   energy__ = Ref{Cdouble}(0.0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_energy(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, energy__::Ref{Cdouble},
                                error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_energy failed with error code", error_code__[])
   end

   return energy__[]
end

function get_forces(gs::GroundStateHandler, label::String)

   forces__ = Matrix{Cdouble}(undef, 3, get_num_atoms(gs))
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_forces(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, forces__::Ref{Cdouble},
                                error_code__::Ref{Cint})::Cvoid 
   if error_code__[] != 0
      error("Sirius.get_forces failed with error code", error_code__[])
   end

   return forces__
end

function get_stress_tensor(gs::GroundStateHandler,  label::String)
   
   stress__ = Matrix{Cdouble}(undef, 3, 3) 
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_stress_tensor(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, stress__::Ref{Cdouble},
                                       error_code__::Ref{Cint})::Cvoid 
   if error_code__[] != 0
      error("Sirius.get_stress_tensor failed with error code", error_code__[])
   end

   return stress__
end

function print_info(ctx::ContextHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_print_info(ctx.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.print_info failed with error code", error_code__[])
   end
end

function print_gs_info(gs::GroundStateHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_print_gs_info(gs.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.print_gs_info failed with error code", error_code__[])
   end
end

end #module
