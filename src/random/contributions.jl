"""
    nullmodel(::Type{C}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, sp::T) where {C <: PermutationConstraint, T}

Returns a probabilistic network where the interaction probabilites for species
`sp` have been replaced by the probabilities given under the null model
specified by the permutation constraint given as its first argument; see
[`nullmodel`](@ref). All interactions that do *not* involve species `sp` have
their probability set to 1.

The original publication [Saavedra2011Strong](@cite) uses [`Degree`] as the
permutation constraint, but this function has been written to be more general.

This network can be passed to [`speciescontribution`](@ref).

###### References

[Saavedra2011Strong](@citet*)
"""
function nullmodel(::Type{C}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, sp::T) where {C <: PermutationConstraint, T}
    R = nullmodel(C, N)
    M = render(Probabilistic{eltype(R.edges)}, N)
    if sp in species(N, 1)
        for suc in species(N, 2)
            M[sp, suc] = R[sp, suc]
        end
    end
    if sp in species(N, 2)
        for pre in species(N, 1)
            M[pre, sp] = R[pre, sp]
        end
    end
    return M
end

"""
TODO
"""
function speciescontribution(::Type{C}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, sp::T, f, args...; replicates=999, kwargs...) where {C <: PermutationConstraint, T}
    R = nullmodel(C, N, sp)
    M = copy(N) # We will update this one in place
    v0 = f(N, args...; kwargs...)
    vx = zeros(typeof(v0), replicates)
    for replicate in 1:replicates
        if sp in species(N, 1)
            i = findfirst(isequal(sp), species(N, 1))
            wv = StatsBase.weights(R[sp,:])
            new_int = StatsBase.sample(axes(N, 2), wv, generality(N)[sp]; replace=false)
            for j in axes(N, 2)
                M[i,j] = j in new_int
            end
        end
        if sp in species(N, 2)
            i = findfirst(isequal(sp), species(N, 2))
            wv = StatsBase.weights(R[:,sp])
            new_int = StatsBase.sample(axes(N, 1), wv, vulnerability(N)[sp]; replace=false)
            for j in axes(N, 1)
                M[j,i] = j in new_int
            end
        end
        vx[replicate] = f(M, args...; kwargs...)
    end
    return (v0 - mean(vx))/std(vx)
end
