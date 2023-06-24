using Revise
using SpeciesInteractionNetworks
import DelimitedFiles

data_path = joinpath(@__DIR__, "data", "nz_stream")

netfiles = filter(endswith(".csv"), readdir(data_path; join=true))

X = DelimitedFiles.readdlm(netfiles[1], ',')
A = permutedims(X[2:end,2:end])
S = [replace(s, r"^(.+) \(N=\d+\)$" => s"\1") for s in X[1,2:end]]
nodes = Unipartite(S)
edges = Binary(Bool.(A))
N = SpeciesInteractionNetwork(nodes, edges)
predecessors(N, "Unidentified detritus")

M = motifs(Unipartite, 3)[1]

sp = rand(species(N))
@profview speciescontribution(Degree, N, sp, complexity)