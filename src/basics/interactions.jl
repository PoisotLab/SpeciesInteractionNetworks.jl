"""
    interactions(N::SpeciesInteractionNetwork)

Returns a vector of all interactions in the network. The interactions are
represented as a a tuple with three elements: the species establishing the
interaction, the species receiving the interaction, and the value of the
interaction.
"""
function interactions(N::SpeciesInteractionNetwork)
    IN = Array{eltype(N)}(undef, length(N))
    i = 1
    for int in N
        IN[i] = int
        i += 1
    end
    return IN
end

@testitem "We can get the interactions in a network" begin
    nodes = Bipartite([:A, :B], [:c, :d])
    edges = Binary([true true; false true])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test interactions(N) == [(:A, :c, true), (:A, :d, true), (:B, :d, true)]
end

@testitem "An empty graph returns and empty interaction list" begin
    nodes = Bipartite([:A, :B], [:c, :d])
    edges = Binary([false false; false false])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test isempty(interactions(N))
end