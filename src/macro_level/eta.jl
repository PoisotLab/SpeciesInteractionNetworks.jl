function η_axis(N::SpeciesInteractionNetwork)
    A = Array(N)
    num = 0.0
    den = 0.0
    for j in 2:richness(N, 1)
        Nj = A[j, :]
        for i in 1:(j - 1)
            Ni = A[i, :]
            num += sum(Ni .* Nj)
            den += min(sum(Ni), sum(Nj))
        end
    end
    return num / den
end

"""
    η(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Probabilistic}})

The η measure of nestedness is a variation of NODF, at the scale of the entire
network. The measure for the entire network is the average of the nestedness of
rows and columns.

###### References

[Bastolla2009architecture](@citet*)
"""
function η(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Probabilistic}})
    return (η(N, 1) + η(N, 2)) / 2.0
end

"""
    η(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Probabilistic}}, dims::Integer)

The η measure of nestedness is a variation of NODF, it can be calculated for
either side of the network (`1` for rows, `2` for columns).

###### References

[Bastolla2009architecture](@citet*)
"""
function η(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Probabilistic}}, dims::Integer)
    dims == 1 && return η_axis(N)
    dims == 2 && return η_axis(permutedims(N))
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

@testitem "η works with known values" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    A = SpeciesInteractionNetwork(nodes, Probabilistic([1.0 0.0 0.0; 0.0 0.1 0.0; 0.0 0.0 0.1]))
    B = SpeciesInteractionNetwork(nodes, Probabilistic([1.0 1.0 1.0; 1.0 0.1 0.0; 1.0 0.0 0.0]))
    C = SpeciesInteractionNetwork(nodes, Probabilistic([1.0 1.0 1.0; 1.0 0.1 0.3; 0.4 0.2 0.0]))
    D = SpeciesInteractionNetwork(nodes, Probabilistic([1.0 1.0 1.0; 0.0 0.1 1.0; 0.0 0.0 1.0]))

    @test η(A) ≈ 0.0
    @test η(B) ≈ 1.0
    @test η(C) ≈ 0.9153846153846155
    @test η(D) ≈ 1.0
end