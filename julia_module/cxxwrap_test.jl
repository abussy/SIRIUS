# Load the module and generate the functions
module SiriusJl
  using CxxWrap
  #TODO: generalize that, this abolute path is very ugly
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

WrappedComm = SiriusJl.CommWrap(comm2f(comm))
SiriusComm = SiriusJl.get_comm(WrappedComm)

@show SiriusJl.initialize(false)

@show ctx = SiriusJl.SimulationContext("./sirius.json", SiriusComm)
@show SiriusJl.initialize(ctx)

k_grid = SiriusJl.get_ngridk(ctx)
k_shift = SiriusJl.get_shiftk(ctx)

@show kps = SiriusJl.KPointSet(ctx, k_grid, k_shift, true)

@show gs = SiriusJl.GroundState(kps)
@show SiriusJl.initial_state(gs)

dtol = SiriusJl.get_density_tol(ctx)
etol = SiriusJl.get_energy_tol(ctx)
itol = SiriusJl.get_initial_tol(ctx)
maxiter = SiriusJl.get_num_dft_iter(ctx)
write_state = true
@show test = SiriusJl.find(gs, dtol, etol, itol, maxiter, write_state)

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

sirius_stress = SiriusJl.stress(gs)
@show smat = SiriusJl.calc_stress_total(sirius_stress)

stress = Matrix{Float64}(undef, 3, 3)
for idx1 = 0:2
   for idx2 = 0:2
      #TODO: hide this in a function somewhere, or even better, in the CxxWrap module
      stress[idx1+1, idx2+1] = SiriusJl.get_element(smat, idx1, idx2)
   end
end
@show stress

@show SiriusJl.finalize(false, true, true)

#Forcing the garbage collector to destroy objects before MPI.Finalize is called (because depends on MPI)
gs = Nothing
kps = Nothing
ctx = Nothing
GC.gc()

MPI.Finalize()
