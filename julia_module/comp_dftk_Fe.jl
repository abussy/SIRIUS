using DFTK
using MKL

@show "Starting calculation ..."

a = 5.42352
lattice = a/2 * [[-1  1  1];
                 [ 1 -1  1];
                 [ 1  1 -1]]
positions = [zeros(3)]

#SIRIUS
Fe = ElementSirius(:Fe; fname=("./Fe_molsim.upf.json"))
atoms     = [Fe]
model = model_SIRIUS(lattice, atoms, positions, ["XC_GGA_X_PBE", "XC_GGA_C_PBE"]; 
                     temperature=0.01, smearing=DFTK.Smearing.FermiDirac(),
                     spin_polarization=:collinear)
UpdateSiriusParams(model, "control", "verbosity", 1)
basis = PlaneWaveBasis(model; Ecut=100, kgrid=[3, 3, 3])
@time begin
   SiriusSCF(basis; density_tol=1.0e-8, energy_tol=1.0e-7, max_niter=5)
   SiriusNlcg(basis)
end
PrintSiriusParams(basis; fname="dftk_sirius.json")
@show GetSiriusEnergy(basis, "total")

#DFTK
Fe = ElementPsp(:Fe; psp=load_psp("./Fe_molsim.upf"))
atoms     = [Fe]
moments = [4]
model = model_PBE(lattice, atoms, positions; temperature=0.01, smearing=DFTK.Smearing.FermiDirac(),
                  spin_polarization=:collinear, magnetic_moments=moments)
@show basis = PlaneWaveBasis(model; Ecut=100, kgrid=[3, 3, 3])
@time begin
   ρ0 = guess_density(basis, moments)
   scfres = self_consistent_field(basis; ρ=ρ0, tol=1e-8, maxiter=100)
end
@show scfres.energies
