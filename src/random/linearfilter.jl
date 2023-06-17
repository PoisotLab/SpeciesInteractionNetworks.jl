"""
    linearfilter(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}; α::Vector{T}=ones(4)) where {T <: AbstractFloat}



###### References

[Stock2017Linear](@citet)
"""
function linearfilter(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}; α::Vector{T}=ones(4)) where {T <: AbstractFloat}
    @assert length(α) == 4
    @assert all(α .>= 0.0 )
    α = α ./ sum(α)

    edges = Probabilistic(
        α[1].*N.edges.edges
            .+ α[2].*mean(N.edges.edges, dims=2)
            .+ α[3].*mean(N.edges.edges, dims=1)
            .+ α[4].*mean(N.edges.edges)
        )
    
    return SpeciesInteractionNetwork(copy(N.nodes), edges)

end

@testitem "The linearfilter function returns the correct type" begin
    nodes = Unipartite([:A, :B, :C, :D, :E, :F])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    R = linearfilter(N)
    @test typeof(R.nodes) <: Unipartite
    @test typeof(R.edges) <: Probabilistic
    @test richness(N) == richness(R)
end

@testitem "The linearfilter with α[4] = 1 is the average only" begin
    import Statistics
    nodes = Unipartite([:A, :B, :C, :D, :E, :F])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    R = linearfilter(N; α=[0.0, 0.0, 0.0, 1.0])
    for interaction in N
        @test R[interaction[1], interaction[2]] == Statistics.mean(N.edges.edges)
    end
end

@testitem "The linearfilter with α[1] = 1 is itself" begin
    import Statistics
    nodes = Unipartite([:A, :B, :C, :D, :E, :F])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    R = linearfilter(N; α=[1.0, 0.0, 0.0, 0.0])
    for interaction in N
        @test R[interaction[1], interaction[2]] == 1.0
    end
end

"""
    nullmodel(::Type{Degree}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

Returns a probabilistic network under the null model that all interactions
occurr with a probability equal to the average of their generality and
vulnerability, *i.e.* proportionally to how many links the species involved are
establishing:
    
``P(i \\rightarrow j) = \\frac{1}{2}\\left(\\frac{\\sum A_{i \\rightarrow \\cdot}}{|B|} + \\frac{\\sum A_{\\cdot \\rightarrow j}}{|T|} \\right)``

###### References

[Bascompte2003nested](@citet)
"""
nullmodel(::Type{Degree}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 1.0, 1.0, 0.0])

"""
    nullmodel(::Type{Connectance}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

Returns a probabilistic network under the null model that all interactions
occurr with a probability equal to the network connectance, *i.e.* in a network
with ``L`` links, and ``|T|`` species on the top and ``|B|`` species on the
bottom,

``P(i \\rightarrow j) = \\frac{L}{|T|\\times |B|}``

Note that for [`Unipartite`](@ref) networks, ``|T| = |B| = |S|`` (the two sides
have the same number of species).

###### References

[Fortuna2006Habitat](@citet)
"""
nullmodel(::Type{Connectance}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 0.0, 0.0, 1.0])

"""
    nullmodel(::Type{Generality}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

Returns a probabilistic network under the null model that interactions happen
proportionally to the generality of the species, *i.e.* their expected number of
outgoing links:

``P(i \\rightarrow j) = \\frac{\\sum A_{i \\rightarrow \\cdot}}{|B|}``

###### References

[Poisot2013Facultative](@citet)

[Weitz2013Phagebacteria](@citet)
"""
nullmodel(::Type{Generality}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 1.0, 0.0, 0.0])

"""
    nullmodel(::Type{Generality}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

Returns a probabilistic network under the null model that interactions happen
proportionally to the vulnerability of the species, *i.e.* their expected number of
incomin links:

``P(i \\rightarrow j) = \\frac{\\sum A_{\\cdot \\rightarrow j}}{|T|}``

###### References

[Poisot2013Facultative](@citet)

[Weitz2013Phagebacteria](@citet)
"""
nullmodel(::Type{Vulnerability}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 0.0, 1.0, 0.0])

@testitem "The degree null model is working" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 1 0; 1 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    R = nullmodel(Degree, N)
    @test typeof(R.nodes) <: Bipartite
    @test typeof(R.edges) <: Probabilistic
    @test R[:A, :a] == 0.5 * (3/3 + 2/3)
    @test R[:B, :b] == 0.5 * (1/3 + 3/3)
    @test R[:C, :c] == 0.5 * (2/3 + 1/3)
end

@testitem "The connectance null model is working" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 1 0; 1 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    R = nullmodel(Connectance, N)
    @test typeof(R.nodes) <: Bipartite
    @test typeof(R.edges) <: Probabilistic
    @test R[:A, :a] == 6 / 9
    @test R[:B, :b] == 6 / 9
    @test R[:C, :c] == 6 / 9
end

@testitem "The generality null model is working" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 1 0; 1 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    R = nullmodel(Generality, N)
    @test typeof(R.nodes) <: Bipartite
    @test typeof(R.edges) <: Probabilistic
    @test R[:A, :a] == 3 / 3
    @test R[:B, :a] == 1 / 3
    @test R[:B, :b] == 1 / 3
    @test R[:B, :c] == 1 / 3
    @test R[:C, :c] == 2 / 3
end

@testitem "The vulnerability null model is working" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 1 0; 1 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    R = nullmodel(Vulnerability, N)
    @test typeof(R.nodes) <: Bipartite
    @test typeof(R.edges) <: Probabilistic
    @test R[:A, :a] == 2 / 3
    @test R[:B, :a] == 2 / 3
    @test R[:B, :b] == 3 / 3
    @test R[:B, :c] == 1 / 3
    @test R[:C, :c] == 1 / 3
end