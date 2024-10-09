"""
    modularity(N::SpeciesInteractionNetwork{<:Partiteness[T], <:Binary}, L::Dict{T,Int}) where {T}

Modularity of a network and its partition. The second argument is a dictionary
where every species of `N` is associated to an integer value representing the
identity of the module. This function returns the same value for bipartite
networks and their unipartite projection. This measure is usually called ``Q``.

Note that in [Newman2006Modularity](@citet), the null model is corrected by
``2m``, where ``m`` is the number of links in the network. This is appropriate
for *undirected* interactions, and so in the code we only divide by ``m``. Note
that this ensure that the modularity matrix has the correct property, *i.e.* its
rows and columns sum to (approx.) 0.

###### References

[Barber2007Modularity](@citet*)

[Newman2006Modularity](@citet*)
"""
function modularity(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Binary}, partition::Dict{T,Int}) where {T}
    
    for s in species(N)
        @assert haskey(partition, s)
    end

    m = links(N)

    ki = [vulnerability(N, s) for s in species(N, 2)]
    ko = [generality(N, s) for s in species(N, 1)]
    si = [partition[s] for s in species(N, 2)]
    so = [partition[s] for s in species(N, 1)]
    P = (ko * ki')/m
    δ = (so .== si')
    B = Array(N) - P
    return sum(B .* δ)/m
end

@testitem "We correctly maximize modularity" begin
    nodes = Bipartite([:A, :B, :C, :D], [:a, :b, :c, :d])
    edges = Binary(zeros(Bool, size(nodes)))
    N = SpeciesInteractionNetwork(nodes, edges)
    N[:A, :a] = true
    N[:A, :b] = true
    N[:B, :a] = true
    N[:B, :b] = true
    N[:C, :c] = true
    N[:C, :d] = true
    N[:D, :d] = true
    N[:D, :c] = true
    partition = Dict([:A => 1, :B => 1, :a => 1, :b => 1, :C => 2, :D => 2, :c => 2, :d => 2])
    @test modularity(N, partition) == 0.5
end
