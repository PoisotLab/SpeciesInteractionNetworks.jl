"""
    degree(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

The degree of a species in a unipartite network is the sum of its generality and
vulnerability. This function *will* count self-links as part of the degree.
"""
function degree(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    deg = Dict([sp => length(predecessors(N, sp))+length(successors(N, sp)) for sp in species(N)])
    for sp in species(N)
        if N[sp,sp]
            deg[sp] -= 1
        end
    end
    return deg
end
