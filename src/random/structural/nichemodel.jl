"""
    structuralmodel(::Type{NicheModel}, species::Integer=10, connectance::AbstractFloat=0.2)

Generate a food web under the niche model with the given number of species, and
an *expected* connectance. Note that by the nature of the niche model algorithm,
only the species richness is guaranteed; there is no guarantee that the
connectance will be correct.

See [`NicheModel`](@ref) for more information about the niche model.
"""
function structuralmodel(::Type{NicheModel}, species::Integer=10, connectance::AbstractFloat=0.2)
    @assert 0.0 < connectance < 0.5
    
    β = 1.0/(2connectance) - 1.0
    
    edges = zeros(Bool, (species, species))
    
    niche = sort(rand(Uniform(0.0, 1.0), species))
    centroids = zeros(species)
    ranges = niche .* rand(Beta(1.0, β), species)

    for species in axes(niche, 1)
        centroids[species] = rand(Uniform(ranges[species]/2, niche[species]))
    end

    for smallest_species in findall(x -> x == minimum(niche), niche)
        ranges[smallest_species] = 0.0
    end

    for consumer in axes(edges,1)
        for resource in axes(edges,2)
            if niche[resource] < centroids[consumer] + 0.5ranges[consumer]
                if niche[resource] > centroids[consumer] - 0.5ranges[consumer]
                    edges[consumer,resource] = true
                end
            end
        end
    end

    edges = Binary(edges)
    nodes = Unipartite(edges)
    return SpeciesInteractionNetwork(nodes, edges)
end

@testitem "We can generate a niche model using structural model" begin
    N = structuralmodel(NicheModel, 12, 0.1)
    @test richness(N) == 12
    @test typeof(N.nodes) <: Unipartite
    @test typeof(N.edges) <: Binary
end

"""
    structuralmodel(::Type{NicheModel}, species::Integer=10, edges::Integer=20)

Generate a food web under the niche model with the given number of species, and
an *expected* numnber of edges. The expected connectance is calculated from the
number of links.

See [`NicheModel`](@ref) for more information about the niche model.
"""
function structuralmodel(::Type{NicheModel}, species::Integer=10, edges::Integer=20)
    @assert species > 0
    @assert (species - 1) < edges < 0.5(species^2)
    return structuralmodel(NicheModel, species, edges/(species^2))
end

@testitem "We can generate a niche model using structural model with a number of links" begin
    N = structuralmodel(NicheModel, 10, 20)
    @test richness(N) == 10
    @test typeof(N.nodes) <: Unipartite
    @test typeof(N.edges) <: Binary
end