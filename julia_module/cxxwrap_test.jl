# Load the module and generate the functions
module SiriusJl
  using CxxWrap
  @wrapmodule(() -> joinpath("/home/bussya/Documents/git/spack/opt/spack/linux-ubuntu23.10-skylake/gcc-13.2.0/sirius-develop-gxaifqlyl44lsm3aylysnyofkahkfp3e/lib/julia_wrapper","libsirius_jl"))

  function __init__()
    @initcxx
  end
end

using MPI
MPI.Init()
comm = MPI.COMM_WORLD

#converting C-style MPI comm to Fortran integer (back and forth, kinda stupid, might want to change that)
comm2f(comm::MPI.Comm) = ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)

@show SiriusJl.initialize(false)

@show ctx = SiriusJl.SimulationContext("./sirius.json")#, comm) looks like issue with type compatibility
@show SiriusJl.initialize(ctx)

#TODO: take these number from the JSON file
k_grid = SiriusJl.R3Vector{Int32}(2, 2, 2)
k_shift = SiriusJl.R3Vector{Int32}(0, 0, 0)

@show kps = SiriusJl.KPointSet(ctx, k_grid, k_shift, true)

@show gs = SiriusJl.GroundState(kps) #TODO: problem of scope => destroyed after MPI.finalize
@show SiriusJl.initial_state(gs)
@show test = SiriusJl.find(gs, 1.0e-6, 1.0e-8, 1.0e-4, 100, true)

energy = SiriusJl.total_energy(gs)
@show energy

sirius_forces = SiriusJl.forces(gs)
@show fmat = SiriusJl.calc_forces_total(sirius_forces)

forces = Matrix{Float64}(undef, 5, 3)

for idx1 = 0:4 #number of atoms, hardcoded for now
   for idx2 = 0:2
      #TODO: hide this in a function somewhere, or even better, in the CxxWrap module
      forces[idx1+1, idx2+1] = SiriusJl.get_element(fmat, 3*idx1 + idx2)
   end
end
@show forces


@show SiriusJl.finalize(false, true, true)

#Forcing the garbage collector to destroy objects before MPI.Finalize is called (because depends on MPI)
gs = Nothing
kps = Nothing
ctx = Nothing
GC.gc()

MPI.Finalize()
