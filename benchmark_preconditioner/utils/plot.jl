using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.YAML
using Base.Meta

results = YAML.load_file("/home/gridsan/hwpang/Jobs/RMS_paper/preconditioner_20230803/results.yml")
results = Dict(mechanism => eval(Meta.parse(result)) for (mechanism, result) in results)

taus = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

numspcs = []
denses = []
preconditionedmaxs = []
preconditionedmins = []
preconditioneddefaults = []
sparsitys = []
besttaus = []
worsttaus = []
for (mech, result) in results
    numspc, numrxn, sparsity = result[1]
    push!(numspcs, numspc)
    push!(sparsitys, sparsity)
    dense = result[2] / 1e9
    push!(denses, dense)
    preconditioned = [result[3][tau] for tau in taus]
    maximum, ind = findmax(preconditioned)
    push!(preconditionedmaxs, maximum / 1e9)
    push!(worsttaus, taus[ind])
    minimum, ind = findmin(preconditioned)
    push!(preconditionedmins, minimum / 1e9)
    push!(besttaus, taus[ind])
    push!(preconditioneddefaults, result[3][1e-3] / 1e9)
end

figure(figsize=(8, 8))
subplot(2, 2, 1)
inds = sortperm(numspcs)
plot(numspcs[inds], denses[inds], "s-", label="Dense", color="C0")
plot(numspcs[inds], preconditionedmaxs[inds], "o--", color="C1", label="Precondtioned - Worst τ")
plot(numspcs[inds], preconditionedmins[inds], "o-", color="C2", label="Preconditioned - Best τ")
plot(numspcs[inds], preconditioneddefaults[inds], "o-.", color="C3", label="Preconditioned - Default τ (1e-3)")
yscale(:log)
xscale(:log)
legend()
xlabel("Number of species")
ylabel("Wall time (s)")

subplot(2, 2, 2)
inds = sortperm(numspcs)
plot(numspcs[inds], besttaus[inds], "o-", color="C2", label="Best")
# plot(numspcs[inds],worsttaus[inds],"o--",color="C1",label="Worst")
yscale(:log)
xscale(:log)
xlabel("Number of species")
ylabel("τ")

subplot(2, 2, 3)
inds = sortperm(sparsitys)
plot(sparsitys[inds] * 100, preconditionedmaxs[inds] ./ denses[inds], "o--", color="C1", label="Worst τ")
plot(sparsitys[inds] * 100, preconditionedmins[inds] ./ denses[inds], "o-", color="C2", label="Best τ")
plot(sparsitys[inds] * 100, preconditioneddefaults[inds] ./ denses[inds], "o-.", color="C3", label="Default τ (1e-3)")
plot(sparsitys[inds] * 100, ones(length(sparsitys)), "k--")
yscale(:log)
xlabel("Sparsity (%)")
ylabel("Wall time ratio (preconditioned/dense)")

subplot(2, 2, 4)
inds = sortperm(sparsitys)
plot(sparsitys[inds] * 100, besttaus[inds], "o-", color="C2", label="Best")
# plot(sparsitys[inds]*100,worsttaus[inds],"d--",color="C2",label="Worst")
yscale(:log)
xlabel("Sparsity (%)")
ylabel("τ")

tight_layout()

savefig("benchmark_preconditioner_taus.pdf")