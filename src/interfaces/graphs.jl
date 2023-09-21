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
See [`reverse`](@ref) for more information.
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

Graphs.has_vertex(N::T, v) where {T <: SpeciesInteractionNetwork} = v in species(N)

@testitem "We can check the existence of a vertex" begin
    import SpeciesInteractionNetworks.Graphs
    edges = Binary(Bool[0 1 0 0; 0 0 1 0; 1 0 0 0; 0 1 1 1])
    nodes = Unipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test Graphs.has_vertex(N, :node_1)
    @test !Graphs.has_vertex(N, :node_1000)
end

function Graphs.has_edge(N::T, s, d) where {T <: SpeciesInteractionNetwork}
    if !(Graphs.has_vertex(N, s) && Graphs.has_vertex(N, d))
        return false
    else
        return !iszero(N[s,d])        
    end
end

@testitem "We can check the existence of an edge" begin
    import SpeciesInteractionNetworks.Graphs
    edges = Binary(Bool[0 1 0 0; 0 0 1 0; 1 0 0 0; 0 1 1 1])
    nodes = Unipartite(edges)
    N = SpeciesInteractionNetwork(nodes, edges)
    @test Graphs.has_edge(N, :node_1, :node_2)
    @test !Graphs.has_edge(N, :node_2, :node_1)
    @test !Graphs.has_edge(N, :node_2000, :node_1)
    @test !Graphs.has_edge(N, :node_2, :node_1000)
    @test !Graphs.has_edge(N, :node_20000, :node_1000)
end

Graphs.inneighbors(N::T, v) where {T <: SpeciesInteractionNetwork} = predecessors(N, v)
Graphs.outneighbors(N::T, v) where {T <: SpeciesInteractionNetwork} = successors(N, v)