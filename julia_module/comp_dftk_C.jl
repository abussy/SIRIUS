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
model = model_SIRIUS(lattice, atoms, positions, ["XC_GGA_X_PBE", "XC_GGA_C_PBE"]; 
                     temperature=0.1, smearing=DFTK.Smearing.FermiDirac(),
                     spin_polarization=:collinear)
UpdateSiriusParams(model, "control", "verbosity", 1)
basis = PlaneWaveBasis(model; Ecut=50, kgrid=[2, 2, 2])
SiriusSCF(basis; density_tol=1.0e-8, energy_tol=1.0e-7, max_niter=4)
SiriusNlcg(basis)
PrintSiriusParams(basis; fname="dftk_sirius.json")
@show GetSiriusEnergy(basis, "total")

#DFTK
C = ElementPsp(:C; psp=load_psp("./C_molsim.upf"))
atoms     = [C, C]
moments = [-1, 1]
model = model_PBE(lattice, atoms, positions; temperature=0.1, smearing=DFTK.Smearing.FermiDirac(),
                  spin_polarization=:collinear, magnetic_moments=moments)
@show basis = PlaneWaveBasis(model; Ecut=50, kgrid=[2, 2, 2])
ρ0 = guess_density(basis, moments)
scfres = self_consistent_field(basis; ρ=ρ0, tol=1e-8, maxiter=100)
@show scfres.energies
