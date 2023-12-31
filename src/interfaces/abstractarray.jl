Base.Array(N::SpeciesInteractionNetwork) = Array(N.edges.edges)
Base.length(N::SpeciesInteractionNetwork) = count(!iszero, N.edges.edges)

@testitem "The length of a network is the number of interactions" begin
    M = rand(Bool, (10, 10))
    N = SpeciesInteractionNetwork{Unipartite, Binary}(M)
    @test length(N) == sum(M)
end

Base.axes(E::Interactions) = axes(E.edges)
Base.axes(E::Interactions, i::Integer) = axes(E.edges, i)

Base.axes(N::SpeciesInteractionNetwork) = axes(N.edges)
Base.axes(N::SpeciesInteractionNetwork, i::Integer) = axes(N.edges, i)

Base.size(E::Interactions) = size(E.edges)
Base.size(E::Interactions, i::Integer) = size(E.edges, i)

Base.size(E::Bipartite) = (length(E.top), length(E.bottom))
Base.size(E::Unipartite) = (length(E.margin), length(E.margin))

Base.size(N::SpeciesInteractionNetwork) = size(N.edges)
Base.size(N::SpeciesInteractionNetwork, i::Integer) = size(N.edges, i)

@testitem "The size of a network is the size of its edges matrix" begin
    M = rand(Bool, (12, 14))
    N = SpeciesInteractionNetwork{Bipartite, Binary}(M)
    @test size(M) == size(N.edges.edges)
    @test size(M, 1) == size(N.edges.edges, 1)
    @test size(M, 2) == size(N.edges.edges, 2)
end

Base.getindex(N::SpeciesInteractionNetwork, args...) = getindex(N.edges, args...)
Base.getindex(E::T, args...) where {T <: Interactions} = getindex(E.edges, args...)
Base.getindex(P::Bipartite{T}, C::CartesianIndex) where {T} = (P.top[C.I[1]], P.bottom[C.I[2]])
Base.getindex(P::Unipartite{T}, C::CartesianIndex) where {T} = (P.margin[C.I[1]], P.margin[C.I[2]])

@testitem "We can get an interaction in the edges of a network by position" begin
    M = rand(Bool, (12, 14))
    E = Binary(M)
    S = Bipartite(E)
    N = SpeciesInteractionNetwork(S, E)
    for i in axes(N, 1)
        for j in axes(E, 2)
            @test E[i, j] == M[i, j]
            @test N[i, j] == M[i, j]
            @test N[i, j] == E[i, j]
        end
    end
end

@testitem "We can get a slice of the network by position" begin
    M = [1 2 3; 4 5 6]
    N = SpeciesInteractionNetwork{Bipartite, Quantitative}(M)
    for i in axes(N, 1)
        @test N[i, :] == M[i, :]
    end
    for j in axes(N, 2)
        @test N[:, j] == M[:, j]
    end
end

function Base.getindex(
    N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions},
    s1::T,
    s2::T,
) where {T}
    i = findfirst(isequal(s1), N.nodes.top)
    j = findfirst(isequal(s2), N.nodes.bottom)
    if isnothing(i)
        throw(ArgumentError("The species $(s1) is not part of the network"))
    end
    if isnothing(j)
        throw(ArgumentError("The species $(s2) is not part of the network"))
    end
    return N[i, j]
end

function Base.getindex(
    N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions},
    s1::T,
    ::Colon,
) where {T}
    i = findfirst(isequal(s1), N.nodes.top)
    if isnothing(i)
        throw(ArgumentError("The species $(s1) is not part of the network"))
    end
    return N[i, :]
end

function Base.getindex(
    N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions},
    ::Colon,
    s2::T,
) where {T}
    j = findfirst(isequal(s2), N.nodes.bottom)
    if isnothing(j)
        throw(ArgumentError("The species $(s2) is not part of the network"))
    end
    return N[:, j]
end

function Base.getindex(
    N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions},
    s1::T,
    s2::T,
) where {T}
    i = findfirst(isequal(s1), N.nodes.margin)
    j = findfirst(isequal(s2), N.nodes.margin)
    if isnothing(i)
        throw(ArgumentError("The species $(s1) is not part of the network"))
    end
    if isnothing(j)
        throw(ArgumentError("The species $(s2) is not part of the network"))
    end
    return N[i, j]
end

function Base.getindex(
    N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions},
    s1::T,
    ::Colon,
) where {T}
    i = findfirst(isequal(s1), N.nodes.margin)
    if isnothing(i)
        throw(ArgumentError("The species $(s1) is not part of the network"))
    end
    return N[i, :]
end

function Base.getindex(
    N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions},
    ::Colon,
    s2::T,
) where {T}
    j = findfirst(isequal(s2), N.nodes.margin)
    if isnothing(j)
        throw(ArgumentError("The species $(s2) is not part of the network"))
    end
    return N[:, j]
end

@testitem "We can index a bipartite network using the species names" begin
    edges = Quantitative([1 2 3; 4 5 6])
    nodes = Bipartite([:A, :B], [:a, :b, :c])
    B = SpeciesInteractionNetwork(nodes, edges)
    @test B[:A, :a] == 1
    @test B[:A, :b] == 2
    @test B[:A, :c] == 3
    @test B[:B, :a] == 4
    @test B[:B, :b] == 5
    @test B[:B, :c] == 6
end

@testitem "We can index a unipartite network using the species names" begin
    edges = Binary([false true; true false])
    nodes = Unipartite([:A, :B])
    U = SpeciesInteractionNetwork(nodes, edges)
    @test U[:A, :] == U[1, :]
    @test U[:, :B] == U[:, 2]
end

@testitem "We can slice a network using species names" begin
    edges = Quantitative([1 2 3; 4 5 6])
    nodes = Bipartite([:A, :B], [:a, :b, :c])
    B = SpeciesInteractionNetwork(nodes, edges)
    @test B[:A, :a] == 1
    @test B[:A, :b] == 2
    @test B[:A, :c] == 3
    @test B[:B, :a] == 4
    @test B[:B, :b] == 5
    @test B[:B, :c] == 6
end

@testitem "We cannot index using a species that is not in the network" begin
    edges = Binary([false true; true false])
    nodes = Unipartite([:A, :B])
    U = SpeciesInteractionNetwork(nodes, edges)
    @test_throws ArgumentError U[:C, :D]
end

function Base.similar(
    N::SpeciesInteractionNetwork{P, E},
) where {P <: Partiteness, E <: Interactions}
    new_edges = E(sparse(zeros(eltype(N.edges), size(N))))
    return SpeciesInteractionNetwork(N.nodes, new_edges)
end

@testitem "We can construct a similar network from a binary network" begin
    N = SpeciesInteractionNetwork{Bipartite, Binary}(rand(Bool, (3, 4)))
    S = similar(N)
    for i in axes(N, 1)
        for j in axes(N, 2)
            @test iszero(S[i, j])
        end
    end
end

@testitem "We can construct a similar network from a quantitative network" begin
    N = SpeciesInteractionNetwork{Unipartite, Quantitative}(rand(Float64, (5, 5)))
    S = similar(N)
    for i in axes(N, 1)
        for j in axes(N, 2)
            @test iszero(S[i, j])
        end
    end
end

@testitem "We can construct a similar network from a probabilistic network" begin
    N = SpeciesInteractionNetwork{Unipartite, Probabilistic}(rand(Float64, (5, 5)))
    S = similar(N)
    for i in axes(N, 1)
        for j in axes(N, 2)
            @test iszero(S[i, j])
        end
    end
end

Base.setindex!(N::SpeciesInteractionNetwork, value, i...) =
    setindex!(N.edges.edges, value, i...)

@testitem "We can use setindex on a network using position indexing" begin
    edges = Quantitative([1 2 3; 4 5 6])
    nodes = Bipartite([:A, :B], [:a, :b, :c])
    B = SpeciesInteractionNetwork(nodes, edges)
    B[1, 1] = 5
    @test B[:A, :a] == 5
    @test B[:A, :b] == 2
end

function Base.setindex!(
    N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions},
    value,
    s1::T,
    s2::T,
) where {T}
    i = findfirst(isequal(s1), N.nodes.top)
    j = findfirst(isequal(s2), N.nodes.bottom)
    if isnothing(i)
        throw(ArgumentError("The species $(s1) is not part of the network"))
    end
    if isnothing(j)
        throw(ArgumentError("The species $(s2) is not part of the network"))
    end
    N[i, j] = value
    return N
end

function Base.setindex!(
    N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions},
    value,
    s1::T,
    s2::T,
) where {T}
    i = findfirst(isequal(s1), N.nodes.margin)
    j = findfirst(isequal(s2), N.nodes.margin)
    if isnothing(i)
        throw(ArgumentError("The species $(s1) is not part of the network"))
    end
    if isnothing(j)
        throw(ArgumentError("The species $(s2) is not part of the network"))
    end
    N[i, j] = value
    return N
end

@testitem "We can use setindex on a network using species indexing" begin
    edges = Quantitative([1 2 3; 4 5 6])
    nodes = Bipartite([:A, :B], [:a, :b, :c])
    B = SpeciesInteractionNetwork(nodes, edges)
    B[:A, :a] = 5
    @test B[:A, :a] == 5
    @test B[:A, :b] == 2
end

@testitem "We can use setindex on a unipartite network using species indexing" begin
    edges = Probabilistic([0.1 0.9; 0.2 0.8])
    nodes = Unipartite([:A, :B])
    B = SpeciesInteractionNetwork(nodes, edges)
    B[:A, :A] = 0.0
    @test iszero(B[:A, :A])
    @test B[:A, :B] == 0.9
end

function Base.getindex(N::SpeciesInteractionNetwork{<:Partiteness{T1}, <:Interactions{T2}}, int::Tuple{T1,T1,T2}) where {T1, T2}
    return N[int[1], int[2]]
end

@testitem "We can access a probabilistic interaction by giving an interaction tuple" begin
    edges = Probabilistic([0.1 0.9; 0.2 0.8])
    nodes = Unipartite([:A, :B])
    B = SpeciesInteractionNetwork(nodes, edges)
    for i in interactions(B)
        @test B[i] == B[i[1], i[2]]
        @test B[i] == i[3]
    end
end

@testitem "We can access a binary interaction by giving an interaction tuple" begin
    edges = Binary(Bool[0 1; 1 1])
    nodes = Unipartite([:A, :B])
    B = SpeciesInteractionNetwork(nodes, edges)
    for i in interactions(B)
        @test B[i] == B[i[1], i[2]]
        @test B[i] == i[3]
    end
end

function Base.setindex!(N::SpeciesInteractionNetwork{<:Partiteness{T1}, <:Interactions{T2}}, val::T2, int::Tuple{T1,T1,T2}) where {T1, T2}
    N[int[1], int[2]] = val
    return N
end

@testitem "We can set an interaction by giving an interaction tuple" begin
    edges = Probabilistic([0.1 0.9; 0.2 0.8])
    nodes = Unipartite([:A, :B])
    B = SpeciesInteractionNetwork(nodes, edges)
    B[(:A, :A, 0.1)] = 1.0
    @test B[:A, :A] == 1.0
end

function Base.permutedims(N::SpeciesInteractionNetwork{<:Bipartite, <:Interactions})
    nodes = Bipartite(N.nodes.bottom, N.nodes.top)
    edges = typeof(N.edges).name.wrapper(permutedims(N.edges.edges))
    return SpeciesInteractionNetwork(nodes, edges)
end

function Base.permutedims(N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions})
    nodes = Unipartite(N.nodes.margin)
    edges = typeof(N.edges).name.wrapper(permutedims(N.edges.edges))
    return SpeciesInteractionNetwork(nodes, edges)
end

@testitem "We can permute the dimensions of a bipartite network" begin
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    edges = Binary(Bool[0 1 1; 1 0 1; 0 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    P = permutedims(N)
    @test sort(species(N)) == sort(species(P))
    @test sort(species(N,1)) == sort(species(P,2))
    @test sort(species(N,2)) == sort(species(P,1))
    for s1 in species(N, 1)
        for s2 in species(N, 2)
            @test N[s1, s2] == P[s2, s1]
        end
    end
end

@testitem "We can permute the dimensions of a unipartite network" begin
    nodes = Unipartite([:A, :B, :C])
    edges = Binary(Bool[0 1 1; 1 0 1; 0 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    P = permutedims(N)
    @test sort(species(N)) == sort(species(P))
    @test sort(species(N,1)) == sort(species(P,2))
    @test sort(species(N,2)) == sort(species(P,1))
    for s1 in species(N, 1)
        for s2 in species(N, 2)
            @test N[s1, s2] == P[s2, s1]
        end
    end
end