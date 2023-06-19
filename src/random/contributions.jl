"""
    nullmodel(::Type{C}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, sp::T) where {C <: PermutationConstraint, T}

Returns a probabilistic network where the interaction probabilites for species
`sp` have been replaced by the probabilities given under the null model
specified by the permutation constraint given as its first argument; see
[`nullmodel`](@ref). All interactions that do *not* involve species `sp` have
their probability set to 1.

The original publication [Saavedra2011Strong](@cite) uses [`Degree`] as the
permutation constraint, but this function has been written to be more general.
It also works on unipartite networks, but this has not been actually *applied*
yet.

This network can be passed to [`speciescontribution`](@ref).

###### References

[Saavedra2011Strong](@citet)
"""
function nullmodel(::Type{C}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, sp::T) where {C <: PermutationConstraint, T}
    R = nullmodel(C, N)
    M = render(Probabilistic{eltype(R.edges)}, N)
    for interaction in N
        if (sp in interaction)
            M[interaction[1], interaction[2]] = R[interaction[1], interaction[2]]
        end
    end
    return M
end

"""
TODO
"""
function speciescontribution(::Type{C}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, sp::T, f, args...; replicates=999, kwargs...) where {C <: PermutationConstraint, T}
    R = nullmodel(C, N, sp)
    v0 = f(N, args...; kwargs...)
    vx = [f(randomdraws(R), args... ;kwargs...) for _ in 1:replicates]
    return (v0 - mean(vx))/std(vx)
end