Base.BroadcastStyle(::Type{T}) where {T <: SpeciesInteractionNetwork} = Broadcast.Style{T}()
Base.BroadcastStyle(::Broadcast.Style{T}, ::S) where {T <: SpeciesInteractionNetwork, S <: Broadcast.BroadcastStyle} = Broadcast.Style{T}()
Base.BroadcastStyle(::S, ::Broadcast.Style{T}) where {T <: SpeciesInteractionNetwork, S <: Broadcast.BroadcastStyle} = Broadcast.Style{T}()

Base.broadcastable(N::T) where {T <: SpeciesInteractionNetwork} = N

find_network(bc::Base.Broadcast.Broadcasted) = find_network(bc.args)
find_network(args::Tuple) = find_network(find_network(args[1]), Base.tail(args))
find_network(x) = x
find_network(::Tuple{}) = nothing
find_network(layer::T, ::Any) where {T <: SpeciesInteractionNetwork} = layer
find_network(::Any, rest) = find_network(rest)

function Base.similar(bc::Base.Broadcast.Broadcasted{Broadcast.Style{T}}, ::Type{ElType}) where {T <: SpeciesInteractionNetwork, ElType}
    network = find_network(bc)
    return similar(network, ElType)
end

function Base.broadcasted(::Broadcast.Style{T}, f, N::T) where {T <: SpeciesInteractionNetwork}
    nodes = copy(N.nodes)
    edges = typeof(N.edges).name.wrapper(f.(Array(N.edges.edges)))
    return SpeciesInteractionNetwork(nodes, edges)
end

function Base.broadcasted(::Broadcast.Style{T}, f, N::T, x) where {T <: SpeciesInteractionNetwork}
    nodes = copy(N.nodes)
    edges = typeof(N.edges).name.wrapper(f.(Array(N.edges.edges), x))
    return SpeciesInteractionNetwork(nodes, edges)
end

function Base.broadcasted(::Broadcast.Style{T}, f, x, N::T) where {T <: SpeciesInteractionNetwork}
    nodes = copy(N.nodes)
    edges = typeof(N.edges).name.wrapper(f.(x, Array(N.edges.edges)))
    return SpeciesInteractionNetwork(nodes, edges)
end

@testitem "We can broadcast over a species interaction network" begin
    M = [1 1 1; 1 1 0; 0 0 1]
    edges = Quantitative(M)
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    Nplus1 = N .+ 1
    for interaction in N
        @test Nplus1[interaction[1], interaction[2]] == interaction[3] + 1.0
    end
    Nplus2 = 2 .+ N
    for interaction in N
        @test Nplus2[interaction[1], interaction[2]] == interaction[3] + 2.0
    end
    NplusM = N .* [2 2 2; 3 3 3; 4 4 4]
    @test NplusM[:A, :A] == 2
    @test NplusM[:A, :B] == 2
    @test NplusM[:A, :C] == 2
    @test NplusM[:B, :A] == 3
    @test NplusM[:B, :B] == 3
    @test NplusM[:B, :C] == 0
    @test NplusM[:C, :A] == 0
    @test NplusM[:C, :B] == 0
    @test NplusM[:C, :C] == 4
end

function Base.broadcasted(
    ::Broadcast.Style{T},
    f,
    U::SpeciesInteractionNetwork{PT, <:Interactions},
    V::SpeciesInteractionNetwork{PT, <:Interactions},
) where {T <: SpeciesInteractionNetwork, PT <: Partiteness}
    @assert sort(species(U)) == sort(species(V))
    nodes = copy(U.nodes)
    ElType = typeof(f(last(first(U)), last(first(V))))
    ItType = Quantitative
    if ElType == Bool
        ItType = Binary
    end
    edges = ItType(zeros(ElType, size(U)))
    N = SpeciesInteractionNetwork(nodes, edges)
    for s1 in species(U)
        for s2 in species(U)
            N[s1, s2] = f(U[s1,s2], V[s1,s2])
        end
    end
    return N
end

@testitem "We can broadcast two networks" begin
    M = [1 1 1; 1 1 0; 0 0 1]
    edges = Quantitative(M)
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    for interaction in interactions(N.+N)
        @test interaction[3] == 2
    end
end

@testitem "We can broadcast unary operations" begin
    M = [1 1 1; 1 1 0; 0 0 1]
    edges = Quantitative(M)
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    for interaction in .-N
        @test interaction[3] == -N[interaction[1], interaction[2]]
    end
end