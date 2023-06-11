function η_axis(N::SpeciesInteractionNetwork)
    S = richness(N, 1)
    n = vec(sum(N.edges.edges, dims=2))
    num = 0.0
    den = 0.0
    for j in 2:S
        Nj = N[j, :]
        for i in 1:(j - 1)
            num += sum(N[i, :] .* Nj)
            den += min(n[i], n[j])
        end
    end
    return num / den
end

"""
    η(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Probabilistic}}, dims::Integer = 0)

The η measure of nestedness is a variation of NODF, it can be calculated at the
scale of the entire network (using `0` as the second argument, which is the
default), or for either side of the network (`1` for rows, `2` for columns).

The measure for the entire network is the average of the nestedness of rows and
columns.

[Bastolla2009Architecture](@cite)
"""
function η(
    N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Probabilistic}},
    dims::Integer = 0,
)
    dims == 1 && return η_axis(N)
    dims == 2 && return η_axis(permutedims(N))
    if iszero(dims)
        return (η(N, 1) + η(N, 2)) / 2.0
    end
    throw(
        ArgumentError(
            "dims can only be 1 (nestedness of rows) or 2 (nestedness of columns), you used $(dims)",
        ),
    )
end

@testitem "The η nestedness of a nested bipartite network is 1" begin
    nodes = Bipartite([:a, :b, :c], [:d, :e, :f])
    edges = Binary(Bool[1 1 1; 1 1 0; 1 0 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test isone(η(N))
    @test isone(η(N,1))
    @test isone(η(N,2))
end