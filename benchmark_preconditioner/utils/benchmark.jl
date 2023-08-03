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

function Reactor_preconditioner(domain::T, y0::Array{T1,1}, tspan::Tuple, interfaces::Z=[]; p::X=DiffEqBase.NullParameters(), forwardsensitivities=false, forwarddiff=false, modelingtoolkit=false, tau=1e-3) where {T<:AbstractDomain,T1<:Real,Z<:AbstractArray,X}

    dydt(dy::X, y::T, p::V, t::Q) where {X,T,Q,V} = dydtreactor!(dy, y, t, domain, interfaces, p=p)
    jacy!(J::Q2, y::T, p::V, t::Q) where {Q2,T,Q,V} = jacobiany!(J, y, p, t, domain, interfaces, nothing)
    jacyforwarddiff!(J::Q2, y::T, p::V, t::Q) where {Q2,T,Q,V} = jacobianyforwarddiff!(J, y, p, t, domain, interfaces, nothing)
    jacp!(J::Q2, y::T, p::V, t::Q) where {Q2,T,Q,V} = jacobianp!(J, y, p, t, domain, interfaces, nothing)
    jacpforwarddiff!(J::Q2, y::T, p::V, t::Q) where {Q2,T,Q,V} = jacobianpforwarddiff!(J, y, p, t, domain, interfaces, nothing)

    @inline function _psetupilu(p::T1, t::T2, u::T3, du::T4, jok::Bool, jcurPtr::T5, gamma::T6, jac!::T7, W::T8, preccache::T9, tau::T10) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
        """
        p: the parameters
        t: the current independent variable
        u: the current state
        du: the current f(u,p,t)
        jok: a bool indicating whether the Jacobian needs to be updated
        jcurPtr: a reference to an Int for whether the Jacobian was updated. jcurPtr[]=true should be set if the Jacobian was updated, and jcurPtr[]=false should be set if the Jacobian was not updated.
        gamma: the gamma of W = M - gamma*J
        """
        if jok

            @. W = 0.0
            jac!(W, u, p, t)
            jcurPtr[] = true

            # W = I - gamma*J
            @. W.nzval = -gamma * W.nzval
            idxs = diagind(W)
            @inbounds @views @. W[idxs] = W[idxs] + 1

            # Build preconditioner on W
            preccache[] = ilu(W, τ=tau)
        end
        nothing
    end
    psetupilu(p::T1, t::T2, u::T3, du::T4, jok::Bool, jcurPtr::T5, gamma::T6) where {T1,T2,T3,T4,T5,T6} = _psetupilu(p, t, u, du, jok, jcurPtr, gamma, jacy!, W::SparseMatrixCSC{Float64,Int64}, preccache::Base.RefValue{IncompleteLU.ILUFactorization{Float64,Int64}}, tau::Float64)
    precilu(z::T1, r::T2, p::T3, t::T4, y::T5, fy::T6, gamma::T7, delta::T8, lr::T9) where {T1,T2,T3,T4,T5,T6,T7,T8,T9} = _precilu(z, r, p, t, y, fy, gamma, delta, lr, preccache)

    J = spzeros(length(y0), length(y0))
    jacy!(J, NaN * ones(length(y0)), p, 0.0)
    @. J.nzval = 1.0

    W = spzeros(length(y0), length(y0))
    jacy!(W, y0, p, 0.0)
    @. W.nzval = -1.0 * W.nzval
    idxs = diagind(W)
    @inbounds @views @. W[idxs] = W[idxs] + 1
    prectmp = ilu(W, τ=tau)
    preccache = Ref(prectmp)

    if (forwardsensitivities || !forwarddiff) && domain isa Union{ConstantTPDomain,ConstantVDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedVDomain,ParametrizedPDomain,ConstantTVDomain,ParametrizedTConstantVDomain,ConstantTAPhiDomain}
        if !forwardsensitivities
            odefcn = ODEFunction(dydt; jac=jacy!, paramjac=jacp!, jac_prototype=J)
        else
            odefcn = ODEFunction(dydt; paramjac=jacp!)
        end
    else
        odefcn = ODEFunction(dydt; jac=jacyforwarddiff!, paramjac=jacpforwarddiff!)
    end
    if forwardsensitivities
        ode = ODEForwardSensitivityProblem(odefcn, y0, tspan, p)
        recsolver = Sundials.CVODE_BDF(linear_solver=:GMRES)
    else
        ode = ODEProblem(odefcn, y0, tspan, p)
        recsolver = Sundials.CVODE_BDF(linear_solver=:GMRES, prec=precilu, psetup=psetupilu, prec_side=1)
    end
    if modelingtoolkit
        sys = modelingtoolkitize(ode)
        jac = eval(ModelingToolkit.generate_jacobian(sys)[2])
        if (forwardsensitivities || !forwarddiff) && domain isa Union{ConstantTPDomain,ConstantVDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedVDomain,ParametrizedPDomain,ConstantTVDomain,ParametrizedTConstantVDomain,ConstantTAPhiDomain}
            odefcn = ODEFunction(dydt; jac=jac, paramjac=jacp!)
        else
            odefcn = ODEFunction(dydt; jac=jac, paramjac=jacpforwarddiff!)
        end
        if forwardsensitivities
            ode = ODEForwardSensitivityProblem(odefcn, y0, tspan, p)
        else
            ode = ODEProblem(odefcn, y0, tspan, p)
        end
    end
    return Reactor(domain, ode, recsolver, forwardsensitivities)
end

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

    b0 = @benchmark sol0 = solve($react0.ode, $Sundials.CVODE_BDF(), abstol=1e-20, reltol=1e-6)
    io = IOBuffer()
    show(io, "text/plain", b0)
    s = String(take!(io))
    println(s)

    bs = Dict{Float64,Float64}()
    @threads for tau in taus

        react = Reactor_preconditioner(domain, y0, (0.0, 1.0); p=p, tau=tau) #Create the reactor object

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