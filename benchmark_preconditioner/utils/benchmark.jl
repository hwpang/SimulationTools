using ReactionMechanismSimulator
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.IncompleteLU
using ReactionMechanismSimulator.LinearAlgebra
using ReactionMechanismSimulator.SparseArrays
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.BenchmarkTools
using Base.Threads

results = Dict()
mechanisms = [
    ("superminimal", "/home/gridsan/hwpang/Software/ReactionMechanismSimulator.jl/src/testing/superminimal.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "H2" => 1.0 / 0.5, "O2" => 1.0, "N2" => 4.0]),),
    ("ethane", "/home/gridsan/hwpang/Software/ReactionMechanismSimulator.jl/src/testing/ethane.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "ethane" => 1.0 / 3.5, "O2" => 1.0, "Ar" => 4.0]),),
    ("sk68", "/home/gridsan/hwpang/Software/SimulationTools/benchmark_preconditioner/mechanism/sk68.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "nc7h16" => 1.0 / 11.0, "n2" => 4.0, "o2" => 1.0]),),
    ("sk88", "/home/gridsan/hwpang/Software/SimulationTools/benchmark_preconditioner/mechanism/sk88.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "nc7h16" => 1.0 / 11.0, "n2" => 4.0, "o2" => 1.0]),),
    ("sk188", "/home/gridsan/hwpang/Software/SimulationTools/benchmark_preconditioner/mechanism/sk188.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "nc7h16" => 1.0 / 11.0, "n2" => 4.0, "o2" => 1.0]),),
    ("PME", "/home/gridsan/hwpang/Software/ReactionMechanismSimulator.jl/src/testing/propyl_methyl_ether.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "PME" => 0.005, "AR" => 0.995]),),
    ("Curran", "/home/gridsan/hwpang/Software/SimulationTools/benchmark_preconditioner/mechanism/Curran.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "nc7h16" => 1.0 / 11.0, "n2" => 4.0, "o2" => 1.0]),),
    ("dodecane1686", "/home/gridsan/hwpang/Software/SimulationTools/benchmark_preconditioner/mechanism/dodecane1686.rms", Dict(["T" => 1000.0, "P" => 10.0e5, "CCCCCCCCCCCC" => 1.0 / 18.5, "N2" => 4.0, "[O][O]" => 1.0]),),
]
taus = [1e-6, 1e-5, 1e-4, 3e-4, 5e-4, 7e-4, 9e-4, 1e-3, 3e-3, 5e-3, 7e-3, 9e-3, 1e-2, 1e-1, 1e0]

function run_test(results, mechanism, taus)
    mechanism_name, input_path, initialconds = mechanism
    J_name = "$(mechanism_name)_J.pdf"

    println(mechanism_name)
    phaseDict = readinput(input_path)
    spcs = phaseDict["phase"]["Species"]
    rxns = phaseDict["phase"]["Reactions"]
    ig = IdealGas(spcs, rxns, name="phase")
    numspcs = length(spcs)
    numrxns = length(rxns)
    println(numspcs)
    println(numrxns)

    initialconds = initialconds
    domain, y0, p = ConstantVDomain(phase=ig, initialconds=initialconds)
    react0 = Reactor(domain, y0, (0.0, 1.0); p=p)

    J = spzeros(length(y0), length(y0))
    react0.ode.f.jac(J, NaN * ones(length(y0)), p, 0.0)
    J.nzval .= 1.0

    sparsity = 1.0 - length(J.nzval) / (size(J)[1] * size(J)[2])

    b0 = @benchmark sol0 = solve($react0.ode, $react0.recommendedsolver, abstol=1e-20, reltol=1e-6)
    io = IOBuffer()
    show(io, "text/plain", b0)
    s = String(take!(io))
    println(s)

    bs = Dict{Float64,Float64}()
    @threads for tau in taus

        react = Reactor_v2(domain, y0, (0.0, 1.0); p=p, tau=tau) #Create the reactor object

        b = @benchmark sol = solve($react.ode, $react.recommendedsolver, abstol=1e-20, reltol=1e-6)
        io = IOBuffer()
        show(io, "text/plain", b)
        s = String(take!(io))

        println(tau)
        println(s)

        bs[tau] = mean(b.times)

    end

    results[mechanism_name] = ((numspcs, numrxns, sparsity), mean(b0.times), bs)
    nothing
end

for mechanism in mechanisms
    run_test(results, mechanism, taus)
end


import YAML
YAML.write_file("benchmark_results.yml", results)