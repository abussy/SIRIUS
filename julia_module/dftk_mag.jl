using Revise
using DFTK                                                                                           
using MKL                                                                                            
using SIRIUS

@show "Starting..."

a = 5.42352  # Bohr
lattice = a / 2 * [[-1  1  1];
                   [ 1 -1  1];
                   [ 1  1 -1]]
atoms     = [ElementPsp(:Fe; psp=load_psp("Fe.upf"; rcut=12.0))]
positions = [zeros(3)]

kgrid = [3, 3, 3] 
Ecut = 32     

magnetic_moments = [4]
#model = model_PBE(lattice, atoms, positions; magnetic_moments, temperature=0.01)
model = model_LDA(lattice, atoms, positions; temperature=0.01)

basis = PlaneWaveBasis(model; Ecut, kgrid)
#ρ0 = guess_density(basis, magnetic_moments)
ρ0 = guess_density(basis)
scfres = self_consistent_field(basis, tol=1e-7; maxiter=100, ρ=ρ0, mixing=SimpleMixing())
scfres.energies

basis = SiriusBasis(model; Ecut, kgrid)
#ρ0 = guess_density(basis)#, magnetic_moments)
scfres = self_consistent_field(basis, tol=1e-7; maxiter=100, ρ=ρ0, mixing=SimpleMixing())
scfres.energies
                                                                                                     
#OLD SILICON TEST
#a = 10.26
#lattice = a / 2 * [[0 1 1.];
#                   [1 0 1.];
#                   [1 1 0.]]
#Si = ElementPsp(:Si; psp=load_psp("Si.upf"; rcut=10.0))
#atoms     = [Si, Si]
#positions = [ones(3)/8, -ones(3)/8]
#magnetic_moments = [2, -2]
#model = model_PBE(lattice, atoms, positions;
#                  magnetic_moments, temperature=0.1, spin_polarization=:collinear)
#
##SIRIUS                                                                                              
#basis = SiriusBasis(model; Ecut=20, kgrid=[2, 2, 2])
#ρ0 = guess_density(basis)
#scfres = self_consistent_field(basis; ρ=ρ0, tol=1.0e-8, maxiter=30, mixing=SimpleMixing())
##@show scfres.eigenvalues[1] .- scfres.eigenvalues[1][1]
##@show scfres.occupation[1]
##@show scfres.energies
#
##DFTK
#println("\n\n\n")
#basis_dftk = PlaneWaveBasis(model; Ecut=20, kgrid=[2, 2, 2])
#ρ0 = guess_density(basis_dftk, magnetic_moments)                                                                   
##ρ0 = guess_density(basis_dftk)                                                                   
#scfres = self_consistent_field(basis_dftk; ρ=ρ0, tol=1e-8, maxiter=30, mixing=SimpleMixing())
##@show scfres.eigenvalues[1] #.- scfres.eigenvalues[1][1]
##@show scfres.occupation[1]
#@show scfres.energies
