function _permutations(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})
    permutations_top = permutations(1:richness(N, 1))
    permutations_bottom = permutations(1:richness(N, 2))
    A = Array(N)
    U = Array{typeof(A),2}(undef, (length(permutations_top), length(permutations_bottom)))
    for (i,pt) in enumerate(permutations_top)
        for (j,pb) in enumerate(permutations_bottom)
            U[i,j] = A[pt,pb]
        end
    end
    return unique(U)
end

function _permutations(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    permutations_sp = permutations(1:richness(N))
    A = Array(N)
    U = Array{typeof(A),1}(undef, length(permutations_sp))
    for (i,ps) in enumerate(permutations_sp)
        U[i] = A[ps,ps]
    end
    return unique(U)
end