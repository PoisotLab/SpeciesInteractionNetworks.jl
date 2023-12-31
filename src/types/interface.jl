Base.eltype(::Bipartite{T}) where {T} = T
Base.eltype(::Unipartite{T}) where {T} = T
Base.eltype(::Binary) = Bool
Base.eltype(::Quantitative{T}) where {T} = T
Base.eltype(::Probabilistic{T}) where {T} = T

@testitem "We can detect the type of a partiteness struct" begin
    part = Bipartite([:A, :B, :C], [:a, :b, :c])
    @test eltype(part) == Symbol

    part = Unipartite(["A", "B", "C"])
    @test eltype(part) == String
end

"""
    species(N::Bipartite)

Returns the list of species in a bipartite list of nodes, as a single vector.
"""
species(N::Bipartite) = vcat(N.top, N.bottom)

"""
    species(N::Unipartite)

Returns the list of species in a unipartite list of nodes, as a single vector.
"""
species(N::Unipartite) = copy(N.margin)

"""
    species(N::SpeciesInteractionNetwork)

Returns the list of species in a network, by calling the `species` method
corresponding to the appropriate species list.
"""
species(N::SpeciesInteractionNetwork) = species(N.nodes)

species(N::Bipartite, dims::Integer) = dims==1 ? copy(N.top) : copy(N.bottom)
species(N::Unipartite, ::Integer) = copy(N.margin)

"""
    species(N::SpeciesInteractionNetwork, dims::Integer)

Returns the list of species on the top (`1` as last argument) or bottom (`2` as
second argument) for the network. For unipartite networks, this will return the
same list of species.
"""
species(N::SpeciesInteractionNetwork, dims::Integer) = species(N.nodes, dims)

@testitem "We can access the species on all sides of the network" begin
    t = [:A, :B, :C]
    b = [:a, :b]

    B = Bipartite(t, b)
    @test species(B, 1) == t
    @test species(B, 2) == b
    @test species(B) == vcat(t,b)
    
    U = Unipartite(t)
    @test species(U, 1) == species(U, 2)
    @test species(U) == species(U, 1)
end

richness(N::Partiteness) = length(species(N))
richness(N::Partiteness, dims::Integer) = length(species(N, dims))

"""
    richness(N::SpeciesInteractionNetwork)

Returns the number of species in a network, measured as the length of the
species items.
"""
richness(N::SpeciesInteractionNetwork) = length(species(N.nodes))

"""
    richness(N::SpeciesInteractionNetwork, dims::Integer)

Returns the number of species in a network, either on the top (`1` as last
argument) or bottom (`2` as last argument), measured as the length of the
species items on this side.
"""
richness(N::SpeciesInteractionNetwork, dims::Integer) = length(species(N.nodes, dims))