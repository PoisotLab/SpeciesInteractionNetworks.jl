"""
    mirror(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

Returns a copy of the network where all interactions are made symmetrical. This
is *not recommended* for daily use, and is intended to use within the
[`spectralradius`](@ref) function.
"""
function mirror(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    M = copy(N)
    for interaction in interactions(N)
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

[Phillips2011structure](@citet*)

[Staniczenko2013ghost](@citet*)
"""
function spectralradius(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary}; correction=:links)
    if iszero(sum(N.edges.edges))
        return NaN
    end
    @assert correction ∈ [:connectance, :links, :none]
    M = mirror(N)
    absolute_real_part = abs.(real.(eigvals(Array(M))))
    L = length(M)
    S = richness(N)
    if correction == :connectance
        return maximum(absolute_real_part)/((L*(S-1))/S)^0.5
    end
    if correction == :links
        return maximum(absolute_real_part)/sqrt((L-sum(diag(Array(N))))/2.0)
    end
    return maximum(absolute_real_part)
end

"""
    spectralradius(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary}; kwargs...)

The bipartite version of [`spectralradius`](@ref) is measured by first
projecting the bipartite network as a unipartite one using [`render`](@ref). The
same options as for the unipartite version are then applied.
"""
function spectralradius(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary}; kwargs...)
    return spectralradius(render(Unipartite, N); kwargs...)
end

@testitem "The values of spectral radius are correct with no correction (test 1)" begin
    edges = Binary(Bool[1 1 1 1; 1 1 1 1; 1 1 1 0; 1 1 0 0; 1 1 0 0; 1 1 0 0])
    nodes = Bipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test spectralradius(N; correction=:none) ≈ 3.82 atol = 0.01
end

@testitem "The values of spectral radius are correct with no correction (test 2)" begin
    edges = Binary(Bool[0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0; 1 1 1 0; 1 0 0 1])
    nodes = Bipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test spectralradius(N; correction=:none) ≈ 3.515 atol = 0.01
end

@testitem "The values of spectral radius are correct with no correction (test 3)" begin
    edges = Binary(Bool[1 1 1 1; 1 1 1 0; 1 1 1 0; 1 1 1 0; 1 1 1 0; 1 0 0 0])
    nodes = Bipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test spectralradius(N; correction=:none) ≈ 3.944 atol = 0.01
end

@testitem "The values of spectral radius are correct with no correction (test 4)" begin
    edges = Binary(Bool[1 0 1 1; 1 1 0 1; 1 1 1 0; 1 1 1 0; 1 1 1 0; 0 1 0 1])
    nodes = Bipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test spectralradius(N; correction=:none) ≈ 3.595 atol = 0.01
end

@testitem "The values of spectral radius are correct with the links corection" begin
    edges = Binary(Bool[1 1 1 1; 1 1 1 0; 1 1 1 0; 1 1 1 0; 1 1 1 0; 1 0 0 0])
    nodes = Bipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test spectralradius(N; correction=:none) ≈ 3.943904 atol = 0.01
    @test spectralradius(N; correction=:links) ≈ 0.956537 atol = 0.01
end
