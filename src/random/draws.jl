"""
    randomdraws(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})

Returns a *binary* network by making random draws from a *probabilistic*
network. Each interaction is considered as an independent Bernoulli trial.

###### References

[Poisot2015structure](@citet)
"""
function randomdraws(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})
    edges = Binary(zeros(Bool, size(N.edges.edges)))
    R = SpeciesInteractionNetwork(copy(N.nodes), edges)
    for interaction in N
        if rand() <= last(interaction)
            R[interaction[1],interaction[2]] = true
        end
    end
    return R
end

@testitem "A random network with all probabilities to one has all interactions" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Probabilistic([1.0 1.0 1.0; 0.0 1.0 0.0; 1.0 0.0 1.0])
    N = SpeciesInteractionNetwork(nodes, edges)
    R = randomdraws(N)
    for interaction in N
        @test R[interaction[1],interaction[2]]
    end
end

@testitem "A random network has only interactions in non-zero positions" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Probabilistic([1.0 0.8 0.1; 0.0 0.2 0.0; 1.0 0.4 1.0])
    N = SpeciesInteractionNetwork(nodes, edges)
    R = randomdraws(N)
    for interaction in R
        @test N[interaction[1],interaction[2]] > 0.0
    end
end