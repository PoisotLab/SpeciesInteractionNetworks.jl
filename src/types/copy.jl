function Base.copy(S::Bipartite)
    return Bipartite(copy(S.top), copy(S.bottom))
end

function Base.copy(S::Unipartite)
    return Unipartite(copy(S.margin))
end

function Base.copy(E::Interactions)
    return (typeof(E).name.wrapper)(copy(E.edges))
end

"""
    Base.copy(N::SpeciesInteractionNetwork)

A copy of a network is created by wrapping together a copy of the nodes and a
copy of the edges.
"""
function Base.copy(N::SpeciesInteractionNetwork)
    return SpeciesInteractionNetwork(copy(N.nodes), copy(N.edges))
end

@testitem "We can make a copy of a network that is not affected by anything" begin
    nodes = Bipartite([:A, :B, :C, :D, :E, :F], [:a, :b, :c, :d, :e, :f, :g, :h])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    P = copy(N)
    for i in 1:10
        swap!(P, Connectance)
    end
    @test P.edges !== N.edges
    @test typeof(P) == typeof(N)
    @test typeof(P.edges) == typeof(N.edges)
end