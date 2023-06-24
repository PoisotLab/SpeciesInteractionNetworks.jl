function findmotif(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary}, M::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    motif_config = _permutations(M)
    sp_combinations = combinations(species(N), richness(M))
    hits = zeros(Bool, length(sp_combinations))
    cursor = 1
    for spcomb in sp_combinations
        hits[cursor] = Array(subgraph(N, spcomb)) in motif_config
        cursor += 1
    end
    return collect(sp_combinations)[hits]
end