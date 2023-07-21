using ReactionMechanismSimulator
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.SparseArrays
using ReactionMechanismSimulator.LinearAlgebra
using ReactionMechanismSimulator.SmoothingSplines
using ReactionMechanismSimulator.StaticArrays
using DataFrames
using CSV

input_file = "chem_liquid_kLA_kH.rms"
Vliq = 1.0
Vgas = 1.0
T = 350 + 273.15
P = 10 * 1e5
C11conc = 0.5 * 740 / (156.31 / 1000)
PDDconc = 0.5 * 856 / (246.43 / 1000)
N2N = P * Vgas / R / T
tf = 30 * 24 * 3600

phaseDict = readinput(input_file)
liqspcs = phaseDict["phase"]["Species"]
liqrxns = phaseDict["phase"]["Reactions"]
solvent = phaseDict["Solvents"][1]
gasspcs = liqspcs
gasrxns = []
liqspcnames = getfield.(liqspcs, :name)
gasspcnames = getfield.(gasspcs, :name)

gas = IdealGas(gasspcs, gasrxns; name="gas");
liq = IdealDiluteSolution(liqspcs, liqrxns, solvent; name="liq", diffusionlimited=true)

initialconds = Dict("C11" => C11conc * Vliq, "PDD" => PDDconc * Vliq, "T" => T, "V" => Vliq)
domainliq, y0liq, pliq = ConstantTVDomain(phase=liq, initialconds=initialconds, constantspecies=["C11", "PDD"]);

initialconds = Dict("N2" => N2N, "T" => T, "P" => P)
domaingas, y0gas, pgas = ConstantTPDomain(phase=gas, initialconds=initialconds);
initialconds = Dict("N2" => N2N, "T" => T, "P" => P)
# inletgas = Inlet(domaingas, initialconds, x -> P * Vgas / R / T)
# outletgas = ConstantVOutlet(domaingas)

masstransferspcnames = getfield.(domainliq.phase.species, :name)
vl, pinter = VaporLiquidMassTransferInternalInterfaceConstantT(domaingas, domainliq, masstransferspcnames);

domains = (domainliq, domaingas)
# interfaces = [vl, inletgas, outletgas]
interfaces = [vl]
react, y0, p = Reactor(domains, (y0liq, y0gas), (0.0, tf), interfaces, (pliq, pgas, pinter));
sol = solve(react.ode, react.recommendedsolver, abstol=1e-18, reltol=1e-6);

df = DataFrame(sol);
CSV.write("Figures/$(save_name).csv", df);

figure()
title("Vapor composition, t=0")
slices = domaingas.indexes[1]:domaingas.indexes[2]
Ns = sol(0)[slices]
N_nz = length(findall(x -> x > 0, Ns))
N = N_nz > 10 ? 10 : N_nz
inds = reverse(sortperm(Ns))[1:N]
pie(Ns[inds], labels=liqspcnames[inds])
savefig("Figures/composition_gas_t0.pdf")

figure()
title("Vapor composition, t=tf")
slices = domaingas.indexes[1]:domaingas.indexes[2]
Ns = sol(tf)[slices]
N_nz = length(findall(x -> x > 0, Ns))
N = N_nz > 10 ? 10 : N_nz
inds = reverse(sortperm(Ns))[1:N]
pie(Ns[inds], labels=liqspcnames[inds])
savefig("Figures/composition_gas_tf.pdf")

figure()
title("Liquid composition, t=0")
slices = domainliq.indexes[1]:domainliq.indexes[2]
Ns = sol(0)[slices]
N_nz = length(findall(x -> x > 0, Ns))
N = N_nz > 10 ? 10 : N_nz
inds = reverse(sortperm(Ns))[1:N]
pie(Ns[inds], labels=liqspcnames[inds])
savefig("Figures/composition_liq_t0.pdf")

figure()
title("Liquid composition, t=tf")
slices = domainliq.indexes[1]:domainliq.indexes[2]
Ns = sol(tf)[slices]
N_nz = length(findall(x -> x > 0, Ns))
N = N_nz > 10 ? 10 : N_nz
inds = reverse(sortperm(Ns))[1:N]
pie(Ns[inds], labels=liqspcnames[inds])
savefig("Figures/composition_liq_tf.pdf")