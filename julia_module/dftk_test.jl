using Revise
using DFTK                                                                                           
using MKL                                                                                            
using SIRIUS

@show "Starting..."
                                                                                                     
a = 3.37                                                                                             
lattice = a * [[0 1 1.];                                                                             
               [1 0 1.];                                                                             
               [1 1 0.]]                                                                             
positions = [ones(3)/8, -ones(3)/8]                                                                  
                                                                                                     
#SIRIUS                                                                                              
C = ElementSirius(:C; fname=("./C_molsim.upf.json"))                                                 
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions, temperature=0.000001)
basis = SiriusBasis(model; Ecut=30, kgrid=[4, 4, 4])                                              

ρ0 = guess_density(basis)
scfres = self_consistent_field(basis; ρ=ρ0, tol=1.0e-10, maxiter=25, mixing=SimpleMixing())
#@show scfres.eigenvalues[1] .- scfres.eigenvalues[1][1]
#@show scfres.occupation[1]
@show scfres.energies.total


#DFTK
println("\n\n\n")
C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"; rcut=10.0))
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions)#, temperature=0.000001)

basis_dftk = PlaneWaveBasis(model; Ecut=30, kgrid=[4, 4, 4])
ρ0 = guess_density(basis_dftk)                                                                   
scfres = self_consistent_field(basis_dftk; ρ=ρ0, tol=1e-10, maxiter=35, mixing=SimpleMixing())
@show scfres.energies.total

# Julia extensions => make SIRIUS an extension to DFTK, such that onlu used if SIRIUS.jl is imported
