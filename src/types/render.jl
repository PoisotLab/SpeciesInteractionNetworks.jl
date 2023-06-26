_nodetype(N::SpeciesInteractionNetwork) = typeof(N.nodes).name.wrapper
_edgetype(N::SpeciesInteractionNetwork) = typeof(N.edges).name.wrapper

"""
    render(::Type{Unipartite}, N::SpeciesInteractionNetwork{<:Bipartite, <:Interactions})

Returns the unipartite projection of a bipartite network. By constructions,
species cannot be shared between levels of bipartite network, so this operation
will always succeed.
"""
function render(::Type{Unipartite}, N::SpeciesInteractionNetwork{<:Bipartite, <:Interactions})
    nodes = Unipartite(species(N,1) ∪ species(N,2))
    edges = _edgetype(N)(zeros(eltype(N.edges), size(nodes)))
    U = SpeciesInteractionNetwork(nodes, edges)
    for sp1 in species(N,1)
        for sp2 in species(N,2)
            if !iszero(N[sp1,sp2])
                U[sp1,sp2] = N[sp1,sp2]
            end
        end
    end
    return U
end

@testitem "We can convert a binary bipartite network to a unipartite network" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 0 1; 1 0 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    M = render(Unipartite, N)
    for interaction in interactions(N)
        @test N[interaction[1], interaction[2]] == M[interaction[1], interaction[2]]
    end
    @test typeof(M.nodes) <: Unipartite
    @test typeof(M.edges) == typeof(N.edges)
end

"""
    render(::Type{Bipartite}, N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions})

Returns the bipartite projection of a unipartite network. By constructions,
species cannot be shared between levels of bipartite network, so this operation
may not always succeed. If it fails, it will throw a
[`BipartiteProjectionFailed`](@ref) exception.
"""
function render(::Type{Bipartite}, N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions})
    top = [s for s in species(N) if isempty(predecessors(N, s))]
    bottom = [s for s in species(N) if isempty(successors(N, s))]
    if !isempty(top ∩ bottom)
        throw(BipartiteProjectionFailed())
    end
    if length(top ∪ bottom) < richness(N)
        throw(BipartiteProjectionFailed())
    end
    nodes = Bipartite(copy(top), copy(bottom))
    edges = _edgetype(N)(zeros(eltype(N.edges), size(nodes)))
    B = SpeciesInteractionNetwork(nodes, edges)
    for interaction in interactions(N)
        B[interaction[1], interaction[2]] = interaction[3]
    end
    return B
end

@testitem "We can convert a binary unipartite network to a bipartite network" begin
    nodes = Unipartite([:A, :B, :a, :b])
    edges = Binary(zeros(Bool, size(nodes)))
    N = SpeciesInteractionNetwork(nodes, edges)
    N[:A, :a] = true
    N[:A, :b] = true
    N[:B, :b] = true
    M = render(Bipartite, N)
    for interaction in interactions(N)
        @test N[interaction[1], interaction[2]] == M[interaction[1], interaction[2]]
    end
    @test typeof(M.nodes) <: Bipartite
    @test richness(M,1) == 2
    @test richness(M,2) == 2
    @test :A ∈ species(M, 1)
    @test :a ∈ species(M, 2)
    @test typeof(M.edges) == typeof(N.edges)
end

@testitem "We cannot convert an overlapping unipartite network to bipartite" begin
    nodes = Unipartite([:A, :B, :a, :b])
    edges = Binary(zeros(Bool, size(nodes)))
    N = SpeciesInteractionNetwork(nodes, edges)
    N[:A, :a] = true
    N[:A, :b] = true
    N[:B, :b] = true
    N[:A, :B] = true
    @test_throws BipartiteProjectionFailed render(Bipartite, N)
end

"""
    render(::Type{Quantitative{T}}, N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions}) where {T <: Number}

Returns a quantitative version of the network, where interaction strengths have
the type `T`. This can be used to convert a quantitative network into a
different number representation.
"""
function render(::Type{Quantitative{T}}, N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions}) where {T <: Number}
    nodes = copy(N.nodes)
    edges = Quantitative(convert.(T, N.edges.edges))
    return SpeciesInteractionNetwork(nodes, edges)
end

@testitem "We can convert a binary network to a quantitative network" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 0 1; 1 0 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    M = render(Quantitative{Float16}, N)
    for interaction in interactions(N)
        @test N[interaction[1], interaction[2]] == M[interaction[1], interaction[2]]
    end
    @test typeof(M.nodes) <: Bipartite
    @test typeof(M.edges) == Quantitative{Float16}
end

@testitem "We can convert a quantitative network to a quantitative network" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Quantitative(Float64[1 2 1; 0 0 1; 1 0 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    M = render(Quantitative{Float32}, N)
    for interaction in interactions(N)
        @test N[interaction[1], interaction[2]] == M[interaction[1], interaction[2]]
    end
    @test typeof(M.nodes) <: Bipartite
    @test typeof(M.edges) == Quantitative{Float32}
end

"""
    render(::Type{Probabilistic{T}}, N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions}) where {T <: AbstractFloat}

Returns a probabilistic version of the network, where interaction probabilities
have the type `T`. This can be used to convert a probabilistic network into a
different number representation.
"""
function render(::Type{Probabilistic{T}}, N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions}) where {T <: AbstractFloat}
    nodes = copy(N.nodes)
    edges = Probabilistic(convert.(T, Array(N)))
    return SpeciesInteractionNetwork(nodes, edges)
end

@testitem "We can turn a binary network into a probabilistic network" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[1 1 1; 0 0 1; 1 0 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    M = render(Probabilistic{Float16}, N)
    for interaction in interactions(N)
        @test isone(M[interaction[1],interaction[2]])
    end
end

"""
    render(::Type{Binary}, N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})

Returns a binary version of the network, where the non-zero interactions are
`true`.
"""
function render(::Type{Binary}, N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})
    nodes = copy(N.nodes)
    edges = Binary(((!) ∘ iszero).(Array(N)))
    return SpeciesInteractionNetwork(nodes, edges)
end