"""
    isdisconnected(N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions}, sp::T) where {T}

Returns `true` if the species has no interaction.
"""
function isdisconnected(N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions}, sp::T) where {T}
    if sp in species(N, 2)
        return isempty(predecessors(N, sp))
    else
        return isempty(successors(N, sp))
    end
end

"""
    isdisconnected(N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions}, sp::T, allow_self_interactions::Bool=true) where {T}

Returns `true` if the species has no interaction; the last argument
(`allow_self_interactions`, defaults to `true`) indicates whether
self-interactions are allowed. If `false`, species that only interact with
themselves are considered to be disconnected.
"""
function isdisconnected(N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions}, sp::T, allow_self_interactions::Bool=true) where {T}
    nei = predecessors(N, sp) ∪ successors(N, sp)
    if !allow_self_interactions
        filter!(!isequal(sp), nei)
    end
    return isempty(nei)
end

"""
    isdegenerate(N::SpeciesInteractionNetwork{<:Bipartite, <:Interactions})

A bipartite network is degenerate if it has species with no interactions.
"""
function isdegenerate(N::SpeciesInteractionNetwork{<:Bipartite, <:Interactions})
    return any([isdisconnected(N, sp) for sp in species(N)])
end

"""
    isdegenerate(N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions}, allow_self_interactions::Bool=true)

A unipartite network is degenerate if it has species with no interactions. In
some cases, it is useful to consider that a species with only self-interactions
is disconnected from the network -- this can be done by using `false` as the
second argument (the default is to *allow* species with only self interactions
to remain).
"""
function isdegenerate(N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions}, allow_self_interactions::Bool=true)
    return any([isdisconnected(N, sp, allow_self_interactions) for sp in species(N)])
end

@testitem "We can identify a network with all disconnected species" begin
    edges = Binary(zeros(Bool, (4, 3)))
    nodes = Bipartite(edges)
    @test isdegenerate(SpeciesInteractionNetwork(nodes, edges))
end

@testitem "We can identify a network with some disconnected species" begin
    edges = Binary(Bool[0 0 0 0; 0 0 0 1; 0 1 1 1; 0 0 0 0])
    nodes = Bipartite(edges)
    @test isdegenerate(SpeciesInteractionNetwork(nodes, edges))
end

@testitem "We can identify a degenerate probabilistic network" begin
    edges = Probabilistic(Float64[0 0 0 0; 0 0 0 0.2; 0 0.9 0.7 0.1; 0 0 0 0])
    nodes = Bipartite(edges)
    @test isdegenerate(SpeciesInteractionNetwork(nodes, edges))
end

@testitem "We can identify a degenerate Quantitative network" begin
    edges = Quantitative(Float64[0 0 0 0; 0 0 0 2; 0 9 7 1; 0 0 0 0])
    nodes = Bipartite(edges)
    @test isdegenerate(SpeciesInteractionNetwork(nodes, edges))
end

@testitem "We can identify a (self-)disconnected species in a unipartite network" begin
    edges = Binary(Bool[1 0 0; 0 0 0; 0 0 1])
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test isdisconnected(N, :A, false)
    @test ~isdisconnected(N, :A, true)
    @test isdisconnected(N, :B, false)
    @test isdisconnected(N, :B, true)
end