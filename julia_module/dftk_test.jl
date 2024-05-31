using Revise
using DFTK                                                                                           
using MKL                                                                                            
using SIRIUS
using LinearAlgebra
                                                                                                     
a = 3.37                                                                                             
lattice = a * [[0 1 1.];                                                                             
               [1 0 1.];                                                                             
               [1 1 0.]]                                                                             
positions = [ones(3)/8, -ones(3)/8]                                                                  
                                                                                                     
#SIRIUS                                                                                              
C = ElementSirius(:C; fname=("./C_molsim.upf.json"))                                                 
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions, temperature=0.000001)
basis = SiriusBasis(model; Ecut=18, kgrid=[2, 2, 2])                                              

ρ0 = guess_density(basis)
scfres = self_consistent_field(basis; ρ=ρ0, tol=1.0e-8, maxiter=15)
@show scfres.ψ[1][1:10, 1]
@show norm(scfres.ψ[1][:, 1])
@show dot(scfres.ψ[1][:, 1], scfres.ψ[1][:, 2])

#energies, H0 = energy_hamiltonian(basis; ρ=ρ0)
#@show SiriusDiagonalize(H0)
#
#function test(num_iter, energies)
#   energies, ham = energy_hamiltonian(basis; ρ=ρ0)
#   SiriusDiagonalize(ham)
#   for i = 2:num_iter
#      ρ = GetSiriusDensity(basis)
#      energies, ham = energy_hamiltonian(basis; ρ=ρ)
#      SiriusDiagonalize(ham)
#      @show energies
#   end
#   return SiriusEnergies(basis)
#end
#
#
#SiriusSCF(basis; density_tol=1.0e-8, energy_tol=1.0e-8, max_niter=100)
#test_energies = test(35, nothing)
#@show test_energies
#@show scfres.energies

#DFTK
println("\n\n\n")
C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"))                                                 
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions, temperature=0.000001)
basis_dftk = PlaneWaveBasis(model; Ecut=18, kgrid=[2, 2, 2])
ρ0 = guess_density(basis_dftk)                                                                   
scfres = self_consistent_field(basis_dftk; ρ=ρ0, tol=1e-8, maxiter=15)
#@show scfres.energies
#@show basis_dftk.kpoints[1].G_vectors[1:10]
@show scfres.ψ[1][1:10, 1]
@show norm(scfres.ψ[1][:, 1])
@show dot(scfres.ψ[1][:, 1], scfres.ψ[1][:, 2])
