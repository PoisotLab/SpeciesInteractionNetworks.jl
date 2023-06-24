"""
    findmotif(M::SpeciesInteractionNetwork{<:Unipartite, <:Binary}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

Search for all groups of species in `N` that match the motif `M` -- note that
the motif to search is the *first* argument.
"""
function findmotif(M::SpeciesInteractionNetwork{<:Unipartite, <:Binary}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    motif_config = _permutations(M)
    A = Array(N)
    sp_combinations = combinations(1:richness(N), richness(M))
    hits = zeros(Bool, length(sp_combinations))
    cursor = 1
    for spcomb in sp_combinations
        hits[cursor] = A[spcomb,spcomb] in motif_config
        cursor += 1
    end
    nads = collect(sp_combinations)[hits]
    T = eltype(N.nodes)
    sps = Vector{Tuple{T,T,T}}(undef, length(nads))
    for (i,nad) in enumerate(nads)
        sps[i] = tuple(species(N)[nad]...)
    end
    return sps
end
