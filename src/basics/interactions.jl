"""
    interactions(N::SpeciesInteractionNetwork)

Returns a vector of all interactions in the network. The interactions are
represented as a a tuple with three elements: the species establishing the
interaction, the species receiving the interaction, and the value of the
interaction.
"""
function interactions(N::SpeciesInteractionNetwork)
    edgelist = Array{eltype(N)}(undef, length(N))
    cursor = 1
    for i in axes(N, 1)
        for j in axes(N, 2)
            if !iszero(N[i,j])
                edgelist[cursor] = (species(N,1)[i],species(N,2)[j],N[i,j]) 
                cursor += 1
            end
        end
    end
    @assert (cursor-1) == length(N) == length(edgelist)
    return edgelist
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