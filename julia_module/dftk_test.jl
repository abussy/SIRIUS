using DFTK                                                                                           
using MKL                                                                                            
                                                                                                     
a = 3.37                                                                                             
lattice = a * [[0 1 1.];                                                                             
               [1 0 1.];                                                                             
               [1 1 0.]]                                                                             
positions = [ones(3)/8, -ones(3)/8]                                                                  
                                                                                                     
#SIRIUS                                                                                              
C = ElementSirius(:C; fname=("./C_molsim.upf.json"))                                                 
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions)                           
basis_sirius = SiriusBasis(model; Ecut=50, kgrid=[2, 2, 2])                                              

ρ0 = guess_density(basis_sirius.PWBasis)
SetSiriusDensity(basis_sirius, ρ0)
rho_backup = ρ0

SiriusSCF(basis_sirius; density_tol=1.0e-8, energy_tol=1.0e-7, max_niter=100)                                 
@show GetSiriusEnergy(basis_sirius, "total")            

#DFTK
C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"))                                                 
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions)                           
basis_dftk = PlaneWaveBasis(model; Ecut=50, kgrid=[2, 2, 2])                                              
ρ0 = guess_density(basis_dftk)                                                                   
#TODO: it seems that rho0 is not necessary, what happens if not provided?
scfres = self_consistent_field(basis_dftk; ρ=ρ0, tol=1e-8, maxiter=100)
@show scfres.energies

scfres = self_consistent_field(basis_dftk; ρ=rho_backup, tol=1e-8, maxiter=100)
@show scfres.energies

#@show rho_backup - ρ0
