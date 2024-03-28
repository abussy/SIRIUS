using DFTK
using MKL

@show "Starting..."

ENV["OMP_NUM_THREADS"] = 4 #TODO: does not seem to work

a = 7.260327248
lattice = a * [[1. 0 0];
               [0 1. 0];
               [0 0 1.]]

Sr = ElementSirius(:Sr; fname="sr_lda_v1.uspp.F.UPF.json")
V  = ElementSirius(:V; fname="v_lda_v1.4.uspp.F.UPF.json")
O  = ElementSirius(:O, fname="o_lda_v1.2.uspp.F.UPF.json")

atoms = [Sr, V, O, O, O]
positions = [[0.45, 0.5, 0.5], [0.0, 0.1, 0.0], [0.5, 0.0, 0.0], [0.0, 0.53, 0.0], [0.0, 0.0, 0.5]]
functionals = ["XC_GGA_X_PBE", "XC_GGA_C_PBE"]

model = model_SIRIUS(lattice, atoms, positions, functionals)
UpdateSiriusParams(model, "control", "verbosity", 1)
UpdateSiriusParams(model, "parameters", "smearing_width", 0.025)
UpdateSiriusParams(model, "parameters", "smearing", "fermi_dirac")
basis = PlaneWaveBasis(model; Ecut=12.0, kgrid=[2, 1, 1])

SiriusSCF(basis)
@show forces = GetSiriusEnergy(basis, "total")
@show forces = GetSiriusForces(basis, "total")
@show stress = GetSiriusStress(basis, "total")

FinalizeSirius(basis)
