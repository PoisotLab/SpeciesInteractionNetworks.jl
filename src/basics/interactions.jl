"""
    interactions(N::SpeciesInteractionNetwork)

Returns a vector of all interactions in the network.

Each interactions is returned an un-named tuple of three elements: the source,
the destination, and the value. For a binary network, for example, an
interaction from `:a` to `:b` will be represented as `(:a,:b,true)`. The type of
the tuple that is returned is given by `eltype(N)`, and is the same as the
output of using iteration on a network.

Note that this method is substantially faster and more memory-efficient than
using iteration (`int for int in N`), for reasons related to the indexing of
non-zero interactions in the underlying data structure of the network.
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