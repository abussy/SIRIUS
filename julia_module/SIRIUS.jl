#This is a tentative attempt at wrapping the SIRIUS C-style API into Julia
module Sirius

export initialize, finalize, create_context, context_handler, free_context_handler

using MPI

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
   error_code__::Int32 = 0
   #TODO: should we look into cfunctions? It looks like most people go that way
   @ccall libpath.sirius_free_object_handler(ctx.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.free_context_handler failed with error code", error_code__)
   end
end

function free_ground_state_handler(gs::ground_state_handler)
   error_code__::Int32 = 0
   @ccall libpath.sirius_free_object_handler(gs.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.free_ground_state_handler failed with error code", error_code__)
   end
end

function free_kpoint_set_handler(kps::kpoint_set_handler)
   error_code__::Int32 = 0
   @ccall libpath.sirius_free_object_handler(kps.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.free_ground_state_handler failed with error code", error_code__)
   end
end

function initialize_context(ctx::context_handler)
   error_code__::Int32 = 0
   @ccall libpath.sirius_initialize_context(ctx.handler_ptr::Ptr{Cvoid}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.initialize_context failed with error code", error_code__)
   end
end


### Simple init/finalize
function initialize(call_mpi_init::Bool)
   call_mpi_init__::Ref{Cuchar} = call_mpi_init
   error_code__::Int32 = 0
   #TODO: is ok to pass error_code__ as an Int32, or should we do Ref{Int32} ? Also prob should use Cint
   @ccall libpath.sirius_initialize(call_mpi_init__::Ref{Cuchar}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.initialize failed with error code", error_code__)
   end
end

function finalize(call_mpi_fin::Bool, call_device_reset::Bool=true, call_fftw_fin::Bool=true)
   call_mpi_fin__::Ref{Cuchar} = call_mpi_fin #TODO: this step is probably not necessary
   call_device_reset__::Ref{Cuchar} = call_device_reset
   call_fftw_fin__::Ref{Cuchar} = call_fftw_fin
   error_code__::Int32 = 0
   @ccall libpath.sirius_finalize(call_mpi_fin__::Ref{Cuchar}, call_device_reset__::Ref{Cuchar}, 
                                  call_fftw_fin__::Ref{Cuchar}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.finalize failed with error code", error_code__)
   end
end

### Utility function that maps an integer to a C MPI_comm
#TODO: -is this complete? 
#      -should we allow C MPI_comms as argument in the C-style API to avoid back and forth?
comm2f(comm::MPI.Comm) = ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)

#TODO: by default, comm_k and comm_band should be the same as comm? (need to verify behavior)
#      according to the Fortran API, if they are not explicitely passed, then they should be null pointers
function create_context(comm::MPI.Comm, ctx::context_handler, comm_k::MPI.Comm, comm_band::MPI.Comm)
   fcomm__::Int32 = comm2f(comm)
   fcomm_k__::Int32 = comm2f(comm_k)
   fcomm_band__::Int32 = comm2f(comm_band)
   error_code__::Int32 = 0
   @ccall libpath.sirius_create_context(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}}, 
                                        fcomm_k__::Int32, fcomm_band__::Int32, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.create_context failed with error code", error_code__)
   end
end

###TESTING a new function that initializes a context from a Json file
function create_context_from_json(comm::MPI.Comm, ctx::context_handler, fname::String)
   fcomm__::Int32 = comm2f(comm)
   error_code__::Int32 = 0
   @ccall libpath.sirius_create_context_from_json(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}},
                                                  fname::Cstring, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.create_context_from_json failed with error code", error_code__)
   end
end

function import_parameters(ctx::context_handler, parameters::String)
   error_code__::Int32 = 0
   @ccall libpath.sirius_import_parameters(ctx.handler_ptr::Ptr{Cvoid}, parameters::Cstring, 
                                           error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.import_parameters failed with error code", error_code__)
   end
end

function set_lattice_vectors(ctx::context_handler, a1::Vector{Float64}, a2::Vector{Float64},
                             a3::Vector{Float64})
   error_code__::Int32 = 0
   @ccall libpath.sirius_set_lattice_vectors(ctx.handler_ptr::Ptr{Cvoid}, a1::Ptr{Cdouble}, a2::Ptr{Cdouble}, 
                                             a3::Ptr{Cdouble}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.set_lattice_vectors failed with error code", error_code__)
   end
end

function add_xc_functional(ctx::context_handler, name::String)
   error_code__::Int32 = 0
   @ccall libpath.sirius_add_xc_functional(ctx.handler_ptr::Ptr{Cvoid}, name::Cstring, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.add_xc_functional failed with error code", error_code__)
   end
end

### Functions for the running of a SCF calculation

function create_kset_from_grid(ctx::context_handler, k_grid::Vector{Int32}, k_shift::Vector{Int32},
                               use_symmetry::Bool, kps::kpoint_set_handler)
   error_code__::Int32 = 0
   use_symmetry__::Ref{Cuchar} = use_symmetry
   @ccall libpath.sirius_create_kset_from_grid(ctx.handler_ptr::Ptr{Cvoid}, k_grid::Ptr{Cint}, 
                                               k_shift::Ptr{Cint}, use_symmetry__::Ref{Cuchar},
                                               kps.handler_ptr::Ptr{Cvoid}, error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.create_kset_from_grid failed with error code", error_code__)
   end
end

function create_ground_state(kps::kpoint_set_handler, gs::ground_state_handler)
   error_code__::Int32 = 0
   @ccall libpath.sirius_create_ground_state(kps.handler_ptr::Ptr{Cvoid}, gs.handler_ptr::Ptr{Cvoid},
                                             error_code__::Int32)::Cvoid
   if error_code__ != 0
      error("Sirius.create_ground_state failed with error code", error_code__)
   end
end


end #module
