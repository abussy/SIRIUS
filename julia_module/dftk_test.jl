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
#positions = [[0.0625, 0.125, 0.125], [-0.0625, -0.125, -0.125],
#             [0.5625, 0.125, 0.125], [0.5-0.0625, -0.125, -0.125]]


#SIRIUS                                                                                              
C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"; rcut=10.0)) #TODO: rcut to be passed to SIRIUS
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions, temperature=0.03)
basis = SiriusBasis(model; Ecut=40, kgrid=[2, 2, 2])                                              

ρ0 = guess_density(basis)
scfres = self_consistent_field(basis; ρ=ρ0, tol=1.0e-8, maxiter=40, mixing=SimpleMixing())
@show scfres.eigenvalues[1]# .- scfres.eigenvalues[1][1]
#@show scfres.occupation[1]
@show scfres.energies

#DFTK
println("\n\n\n")
C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"; rcut=10.0))
atoms     = [C, C]                                                                                   
model = model_PBE(lattice, atoms, positions, temperature=0.03)

basis_dftk = PlaneWaveBasis(model; Ecut=40, kgrid=[2, 2, 2])
ρ0 = guess_density(basis_dftk)                                                                   
scfres = self_consistent_field(basis_dftk; ρ=ρ0, tol=1e-8, maxiter=40, mixing=SimpleMixing())
@show scfres.eigenvalues[1] #.- scfres.eigenvalues[1][1]
#@show scfres.occupation[1]
@show scfres.energies
