using Revise
using SpeciesInteractionNetworks
import DelimitedFiles
import Statistics

data_path = joinpath(@__DIR__, "data", "nz_stream")

netfiles = filter(endswith(".csv"), readdir(data_path; join=true))

for nf in netfiles
    X = DelimitedFiles.readdlm(nf, ',')
    A = permutedims(X[2:end,2:end])
    S = [replace(s, r"^(.+) \(N=\d+\)$" => s"\1") for s in X[1,2:end]]
    nodes = Unipartite(S)
    edges = Binary(Bool.(A))
    N = SpeciesInteractionNetwork(nodes, edges)
    @info last(splitpath(nf))
    R = [structuralmodel(NicheModel, richness(N), connectance(N)) for i in 1:150]
    r = complexity.(R)
    m, v = Statistics.mean(r), Statistics.std(r)
    x = complexity(N)
    @info " -> $((x - m)/v)"
end

predecessors(N, "Unidentified detritus")

M = motifs(Unipartite, 3)[1]

sp = rand(species(N))
@profview speciescontribution(Degree, N, sp, complexity)