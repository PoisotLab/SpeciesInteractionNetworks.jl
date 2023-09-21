"""
    Base.reverse(N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})

Returns a copy of the network, in which the interactions have been flipped. In
other words, an interaction A → B is now B → A. This maintains the nature of the
interaction, and works for loops, self-edges etc.
"""
function Base.reverse(N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})
    M = copy(N)
    reverse!(M)
end

"""
    Base.reverse!(N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})

Modifies the network given as its argument so that the interactions are flipped.
See [`reverese`](@ref) for more information.
"""
function Base.reverse!(N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})
    for interaction in interactions(N)
        N[interaction[2], interaction[1]], N[interaction[1], interaction[2]] = N[interaction[1], interaction[2]], N[interaction[2], interaction[1]]
    end
    return N
end

@testitem "We can reverse a network" begin
    edges = Binary(Bool[0 1 0 0; 0 0 1 0; 1 0 0 0; 0 1 1 1])
    nodes = Unipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    R = reverse(N)
    for interaction in interactions(R)
        @test N[interaction[2], interaction[1]] == interaction[3]
    end
end