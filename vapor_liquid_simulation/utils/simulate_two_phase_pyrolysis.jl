using ReactionMechanismSimulator
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.SparseArrays
using ReactionMechanismSimulator.LinearAlgebra
using ReactionMechanismSimulator.SmoothingSplines
using ReactionMechanismSimulator.StaticArrays
using DataFrames
using CSV

input_file = "/home/gridsan/hwpang/Software/SimulationTools/vapor_liquid_simulation/mechanism/chem_solvation_kLA_kH.rms"
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
inletgas = Inlet(domaingas, initialconds, x -> P * Vgas / R / T)
outletgas = MaintainingVOutlet(domaingas)

masstransferspcnames = getfield.(domainliq.phase.species, :name)
vl, pinter = VaporLiquidMassTransferInternalInterfaceConstantT(domaingas, domainliq, masstransferspcnames);

domains = (domainliq, domaingas)
interfaces = [vl, inletgas, outletgas]
react, y0, p = Reactor(domains, (y0liq, y0gas), (0.0, tf), interfaces, (pliq, pgas, pinter));
sol = solve(react.ode, react.recommendedsolver, abstol=1e-18, reltol=1e-6);

if !isdir("Figures")
    mkdir("Figures")
end

df = DataFrame(sol);
CSV.write("Figures/$(save_name).csv", df);


ts = range(0.0, tf, length=50)

figure(figsize=(12, 4))
subplot(1, 2, 1)
title("Vapor")
plotmolefractions(ssys.sims[2], ssys.sol.t[end-1], tol=7e-2)

subplot(1, 2, 2)
title("Liquid")
plotmolefractions(ssys.sims[1], ssys.sol.t[end-1], tol=7e-2)
savefig("./Figures/two_phase.pdf", bbox_inches="tight")

# import ReactionMechanismSimulator.getfluxdiagram
# import ReactionMechanismSimulator.makefluxdiagrams
# using ReactionMechanismSimulator: getdomainsize, getcolor
# using ReactionMechanismSimulator.PyCall

# @pyimport pydot
# function getfluxdiagram(bsol, t; centralspecieslist=Array{String,1}(), superimpose=false,
#     maximumnodecount=50, maximumedgecount=50, concentrationtol=1e-6, speciesratetolerance=1e-6,
#     maximumnodepenwidth=10.0, maximumedgepenwidth=10.0, radius=1, centralreactioncount=-1, outputdirectory="fluxdiagrams",
#     colorscheme="viridis", removeunconnectednodes=true)

#     fd = makefluxdiagrams(bsol, [t]; centralspecieslist=centralspecieslist, superimpose=superimpose,
#         maximumnodecount=maximumnodecount, maximumedgecount=maximumedgecount, concentrationtol=concentrationtol,
#         speciesratetolerance=speciesratetolerance, maximumnodepenwidth=maximumnodepenwidth,
#         maximumedgepenwidth=maximumedgepenwidth, radius=radius, centralreactioncount=centralreactioncount,
#         outputdirectory=outputdirectory, colorscheme=colorscheme, removeunconnectednodes=removeunconnectednodes)

#     return getdiagram(fd, 1)
# end

# function makefluxdiagrams(bsol, ts; centralspecieslist=Array{String,1}(), superimpose=false,
#     maximumnodecount=50, maximumedgecount=50, concentrationtol=1e-6, speciesratetolerance=1e-6,
#     maximumnodepenwidth=10.0, maximumedgepenwidth=10.0, radius=1, centralreactioncount=-1, outputdirectory="fluxdiagrams",
#     colorscheme="viridis", removeunconnectednodes=false)

#     if hasproperty(bsol, :domain)
#         specieslist = bsol.domain.phase.species
#         reactionlist = bsol.domain.phase.reactions
#     else
#         specieslist = bsol.species
#         reactionlist = bsol.reactions
#     end

#     speciesnamelist = getfield.(specieslist, :name)
#     numspecies = length(specieslist)

#     if !isdir(outputdirectory)
#         mkdir(outputdirectory)
#     end

#     concs = reduce(vcat, [concentrations(bsol, t) for t in ts])

#     reactionrates = reduce(vcat, [rates(bsol, t) for t in ts])

#     drawspecies(specieslist)
#     speciesdirectory = joinpath(pwd(), "species")

#     #find central species
#     centralspeciesindices = Array{Int64,1}()
#     if length(centralspecieslist) != 0
#         for (j, centralspecies) in enumerate(centralspecieslist)
#             for (i, species) in enumerate(specieslist)
#                 if species.name == centralspecies
#                     push!(centralspeciesindices, i)
#                     break
#                 end
#             end
#             if length(centralspeciesindices) != j
#                 throw(error("Central species $centralspecies could not be found"))
#             end
#         end
#     end

#     speciesrates = zeros(numspecies, numspecies, length(ts))
#     for (index, reaction) in enumerate(reactionlist)
#         rate = reactionrates[index, :]
#         if length(reaction.pairs[1]) > 1
#             pairs = reaction.pairs
#         else
#             pairs = getpairs(reaction)
#         end
#         for (reactant, product) in pairs
#             reactantindex = findfirst(y -> y == reactant, speciesnamelist)
#             productindex = findfirst(y -> y == product, speciesnamelist)
#             speciesrates[reactantindex, productindex, :] .+= rate
#             speciesrates[productindex, reactantindex, :] .-= rate
#         end
#     end

#     maxconcs = maximum(concs, dims=2)
#     maxconcentration = maximum(maxconcs)

#     maxreactionrates = maximum(abs.(reactionrates), dims=2)

#     maxspeciesrates = maximum(abs.(speciesrates), dims=3)
#     maxspeciesrate = maximum(maxspeciesrates)

#     speciesindex = sortperm(reshape(maxspeciesrates, numspecies^2), rev=true)

#     nodes = Array{Int64,1}()
#     edges = []
#     if !superimpose && length(centralspecieslist) != 0
#         for centralspeciesindex in centralspeciesindices
#             push!(nodes, centralspeciesindex)

#             if radius == 0
#                 for reaction in reactionlist
#                     if length(reaction.pairs[1]) > 1
#                         pairs = reaction.pairs
#                     else
#                         pairs = getpairs(reaction)
#                     end
#                     for (reactant, product) in pairs
#                         rindex = findfirst(y -> y == reactant, speciesnamelist)
#                         pindex = findfirst(y -> y == product, speciesnamelist)
#                         if rindex in nodes && pindex in nodes
#                             if !((rindex, pindex) in edges) && !((pindex, rindex) in edges)
#                                 push!(edges, (rindex, pindex))
#                             end
#                         end
#                     end
#                 end
#             else
#                 addadjacentnodes!(centralspeciesindex, nodes, edges, reactionlist,
#                     maxreactionrates, maxspeciesrates, centralreactioncount, radius, Array{Int64,1}(), speciesnamelist)
#             end
#         end
#     else
#         for i = 1:numspecies^2
#             productindex = div(speciesindex[i], numspecies) + 1
#             reactantindex = rem(speciesindex[i], numspecies)
#             if reactantindex == 0
#                 reactantindex = numspecies
#                 productindex -= 1
#             end
#             if reactantindex > productindex
#                 continue
#             end
#             if maxspeciesrates[reactantindex, productindex] == 0
#                 break
#             end
#             if !(reactantindex in nodes) && length(nodes) < maximumnodecount
#                 push!(nodes, reactantindex)
#             end
#             if !(productindex in nodes) && length(nodes) < maximumnodecount
#                 push!(nodes, productindex)
#             end
#             if !((reactantindex, productindex) in edges) && !((productindex, reactantindex) in edges)
#                 push!(edges, (reactantindex, productindex))
#             end
#             if length(nodes) > maximumnodecount
#                 break
#             end
#             if length(edges) >= maximumedgecount
#                 break
#             end
#         end
#         if superimpose && length(centralspecieslist) > 0
#             nodescopy = nodes[:]
#             for centralspeciesindex in centralspeciesindices
#                 if !(centralspeciesindex in nodes)
#                     push!(nodes, centralspeciesindex)
#                     addadjacentnodes!(centralspeciesindex, nodes, edges, reactionlist,
#                         maxreactionrates, maxspeciesrates, centralreactioncount, -1, nodescopy, speciesnamelist)
#                 end
#             end
#         end
#     end



#     if removeunconnectednodes
#         if length(ts) > 1
#             error("cannot use removeunconnectednodes for length(ts)>1")
#         end
#         minspeciesrate = Inf
#         maxspcrate = -Inf
#         for index in 1:length(edges)
#             reactantindex, productindex = edges[index]
#             sprate = abs(speciesrates[reactantindex, productindex, 1])
#             if sprate > maxspcrate
#                 maxspcrate = sprate
#             end
#         end
#         connectedinds = []
#         for edge in edges
#             reactantindex, productindex = edge
#             speciesrate = speciesrates[reactantindex, productindex, 1] / maxspcrate
#             for ind in edge
#                 if abs(speciesrate) > speciesratetolerance && !(ind in connectedinds)
#                     push!(connectedinds, ind)
#                 end
#             end
#         end
#         filter!(x -> (x in connectedinds), nodes)
#     end

#     graph = pydot.Dot("flux_diagram", graph_type="digraph", overlap="false")
#     graph.set_rankdir("LR")
#     graph.set_fontname("sans")
#     graph.set_fontsize("10")

#     for index in nodes
#         species = specieslist[index]
#         node = pydot.Node(name=species.name)
#         node.set_penwidth(maximumnodepenwidth)
#         graph.add_node(node)

#         speciesindex = string(species.name, ".png")
#         imagepath = ""
#         if !isdir(speciesdirectory)
#             continue
#         end
#         for (root, dirs, files) in walkdir(speciesdirectory)
#             for f in files
#                 if f == speciesindex
#                     imagepath = joinpath(root, f)
#                     break
#                 end
#             end
#         end
#         if isfile(imagepath)
#             node.set_image(imagepath)
#             node.set_label(" ")
#         end
#     end

#     for (reactantindex, productindex) in edges
#         if reactantindex in nodes && productindex in nodes
#             reactant = specieslist[reactantindex]
#             product = specieslist[productindex]
#             edge = pydot.Edge(reactant.name, product.name)
#             edge.set_penwidth(maximumedgepenwidth)
#             graph.add_edge(edge)
#         end
#     end

#     graph = pydot.graph_from_dot_data(graph.create_dot(prog="dot"))[1]

#     for t in 1:length(ts)
#         slope = -maximumnodepenwidth / log10(concentrationtol)
#         for index in nodes
#             species = specieslist[index]
#             if occursin(r"^[a-zA-Z0-9_]*$", species.name)
#                 species_string = species.name
#             else
#                 species_string = string("\"", species.name, "\"")
#             end

#             node = graph.get_node(species_string)[1]
#             concentration = concs[index, t] / maxconcentration
#             if concentration < concentrationtol
#                 penwidth = 0.0
#             else
#                 penwidth = round((slope * log10(concentration) + maximumnodepenwidth) * 1.0e3) / 1.0e3
#             end

#             node.set_penwidth(penwidth)
#         end

#         slope = -maximumedgepenwidth / log10(speciesratetolerance)
#         maxspcrate = -Inf
#         minspeciesrate = Inf
#         for index in 1:length(edges)
#             reactantindex, productindex = edges[index]
#             sprate = abs(speciesrates[reactantindex, productindex, t])
#             if minspeciesrate > sprate
#                 minspeciesrate = sprate
#             end
#             if sprate > maxspcrate
#                 maxspcrate = sprate
#             end
#         end
#         minspeciesrate = max(minspeciesrate, maxspcrate * speciesratetolerance)

#         for index in 1:length(edges)
#             reactantindex, productindex = edges[index]
#             if reactantindex in nodes && productindex in nodes
#                 reactant = specieslist[reactantindex]
#                 product = specieslist[productindex]

#                 if occursin(r"^[a-zA-Z0-9_]*$", reactant.name)
#                     reactant_string = reactant.name
#                 else
#                     reactant_string = string("\"", reactant.name, "\"")
#                 end

#                 if occursin(r"^[a-zA-Z0-9_]*$", product.name)
#                     product_string = product.name
#                 else
#                     product_string = string("\"", product.name, "\"")
#                 end

#                 edge = graph.get_edge(reactant_string, product_string)[1]

#                 speciesrate = speciesrates[reactantindex, productindex, t] / maxspcrate
#                 if speciesrate < 0
#                     edge.set_dir("back")
#                     speciesrate = -speciesrate
#                 else
#                     edge.set_dir("forward")
#                 end

#                 if speciesrate < speciesratetolerance
#                     penwidth = 0.0
#                     edge.set_dir("none")
#                 else
#                     penwidth = round((slope * log10(speciesrate) + maximumedgepenwidth) * 1.0e3) / 1.0e3
#                 end

#                 edge.set_penwidth(penwidth)
#                 edge.set_color(getcolor(speciesrates[reactantindex, productindex, t], maxspcrate, minspeciesrate, colorscheme))

#             end
#         end

#         if ts[t] == 0.0
#             label = "t = 0 s"
#         else
#             tval = log10(ts[t])
#             label = "t = 10^$tval s"
#         end

#         graph.set_label(label)
#         graph.write_dot(joinpath(outputdirectory, "flux_diagram_$t.dot"))
#         graph.write_png(joinpath(outputdirectory, "flux_diagram_$t.png"))
#         graph.write_svg(joinpath(outputdirectory, "flux_diagram_$t.svg"))
#         graph.write_pdf(joinpath(outputdirectory, "flux_diagram_$t.pdf"))
#     end
#     return FluxDiagram(ts, outputdirectory)
# end
# function rates(ssys::Q, t::X) where {Q<:SystemSimulation,X<:Real}
#     rts = zeros(length(ssys.reactions))
#     domains = getfield.(ssys.sims, :domain)
#     Nrxns = sum([length(sim.domain.phase.reactions) for sim in ssys.sims]) + sum([hasproperty(inter, :reactions) ? length(inter.reactions) : 0 for inter in ssys.interfaces])
#     Nspcs = sum([length(sim.domain.phase.species) for sim in ssys.sims])
#     cstot = zeros(Nspcs)
#     vns = Array{Any,1}(undef, length(domains))
#     vcs = Array{Any,1}(undef, length(domains))
#     vT = Array{Any,1}(undef, length(domains))
#     vP = Array{Any,1}(undef, length(domains))
#     vV = Array{Any,1}(undef, length(domains))
#     vC = Array{Any,1}(undef, length(domains))
#     vN = Array{Any,1}(undef, length(domains))
#     vmu = Array{Any,1}(undef, length(domains))
#     vkfs = Array{Any,1}(undef, length(domains))
#     vkrevs = Array{Any,1}(undef, length(domains))
#     vHs = Array{Any,1}(undef, length(domains))
#     vUs = Array{Any,1}(undef, length(domains))
#     vGs = Array{Any,1}(undef, length(domains))
#     vdiffs = Array{Any,1}(undef, length(domains))
#     vCvave = Array{Any,1}(undef, length(domains))
#     vphi = Array{Any,1}(undef, length(domains))
#     index = 1
#     for (k, sim) in enumerate(ssys.sims)
#         vns[k], vcs[k], vT[k], vP[k], vV[k], vC[k], vN[k], vmu[k], vkfs[k], vkrevs[k], vHs[k], vUs[k], vGs[k], vdiffs[k], vCvave[k], vphi[k] = calcthermo(sim.domain, ssys.sol(t), t)
#         cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] = vcs[k]
#         rts[index:index+length(vkfs[k])-1] .= getrates(sim.domain.rxnarray, cstot, vkfs[k], vkrevs[k]) .* getdomainsize(sim, t)
#         index += length(vkfs[k])
#     end
#     for inter in ssys.interfaces
#         if hasproperty(inter, :reactions)
#             kfs, krevs = getkfskrevs(inter, vT[inter.domaininds[1]], vT[inter.domaininds[2]], vphi[inter.domaininds[1]], vphi[inter.domaininds[2]], vGs[inter.domaininds[1]], vGs[inter.domaininds[2]], cstot)
#             rts[index:index+length(kfs)-1] = getrates(inter.rxnarray, cstot, kfs, krevs) .* inter.A
#             index += length(kfs)
#         end
#     end
#     return rts
# end
getfluxdiagram(ssys, tf, concentrationtol=1e-3)