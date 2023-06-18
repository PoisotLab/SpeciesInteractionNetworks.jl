"""
    mirror(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

Returns a copy of the network where all interactions are made symmetrical. This
is *not recommended* for daily use, and is intended to use within the
[`spectralradius`](@ref) function.
"""
function mirror(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    M = copy(N)
    for interaction in N
        M[interaction[2], interaction[1]] = true
    end
    return M
end

"""
    spectralradius(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary}; correction=:links)

The spectral radius of a unipartite is a conceptual equivalent to nestedness
[Staniczenko2013ghost](@cite). It is defined as the absolute value of the
largest real part of the eigenvalues of the *undirected* adjacency matrix.

There are a number of corrections available through the `correction` keyword.

The default correction is `:links` as in [Staniczenko2013ghost](@citet), where
the values are divided by the square root of the number of links, *excluding the
self-interactions*.

Using the `:connectance` correction follows the version of
[Phillips2011structure](@citet), where the values are divided by the square root
of ``(L\\times(S-1))S^{-1}`` (this is not *quite* connectance, but the point is
that this version is corrected for network size and order).

Using `:none` returns the raw values, and for the sake of comparisons across
networks, it is advised not to use it. It is included for cases where the
networks to compare have the same number of species and interactions, as in this
case it is appropriate and slightly faster than other corrections.

###### References

[Phillips2011structure](@citet)

[Staniczenko2013ghost](@citet)
"""
function spectralradius(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary}; correction=:links)
    if iszero(sum(N.edges.edges))
        return NaN
    end
    @assert correction ∈ [:connectance, :links, :none]
    M = mirror(N)
    absolute_real_part = abs.(real.(eigvals(Array(M.edges.edges))))
    if correction == :connectance
        return maximum(absolute_real_part)/((length(N)*(richness(N)-1))/richness(N))^0.5
    end
    if correction == :links
        return maximum(absolute_real_part)/sqrt((length(N)-sum(diag(Array(N.edges.edges))))/2.0)
    end
    return maximum(absolute_real_part)
end