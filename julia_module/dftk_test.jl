using DFTK                                                                                           
using MKL                                                                                            
using SIRIUS
                                                                                                     
a = 3.37                                                                                             
lattice = a * [[0 1 1.];                                                                             
               [1 0 1.];                                                                             
               [1 1 0.]]                                                                             
positions = [ones(3)/8, -ones(3)/8]                                                                  
                                                                                                     
#SIRIUS                                                                                              
C = ElementSirius(:C; fname=("./C_molsim.upf.json"))                                                 
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions, temperature=0.1)                           
basis_sirius = SiriusBasis(model; Ecut=50, kgrid=[2, 2, 2])                                              

ρ0 = guess_density(basis_sirius.PWBasis)
#H0 = SiriusHamiltonian(basis_sirius, ρ0)
#energies = SiriusEnergies(basis_sirius)

function test(num_iter, energies)
   energies, ham = energy_hamiltonian(basis_sirius; ρ=ρ0)
   SiriusDiagonalize(ham)
   for i = 2:num_iter
      ρ = GetSiriusDensity(basis_sirius)
      energies, ham = energy_hamiltonian(basis_sirius; ρ=ρ)
      SiriusDiagonalize(ham)
   end
   return SiriusEnergies(basis_sirius)
end


SiriusSCF(basis_sirius; density_tol=1.0e-8, energy_tol=1.0e-8, max_niter=100)
test_energies = test(35, nothing)

#DFTK
#C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"))                                                 
#atoms     = [C, C]                                                                                   
#model = model_PBE(lattice, atoms, positions) #, temperature=0.1)
#basis_dftk = PlaneWaveBasis(model; Ecut=50, kgrid=[2, 2, 2])                                              
#ρ0 = guess_density(basis_dftk)                                                                   
##TODO: it seems that rho0 is not necessary, what happens if not provided?
#scfres = self_consistent_field(basis_dftk; ρ=ρ0, tol=1e-8, maxiter=100)
#@show scfres.energies
#
#test_energies
