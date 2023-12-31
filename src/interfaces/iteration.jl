function Base.IteratorSize(::Type{T}) where {T<:SpeciesInteractionNetwork}
    return Base.HasLength()
end

Base.eltype(S::SpeciesInteractionNetwork) = Tuple{eltype(S.nodes), eltype(S.nodes), eltype(S.edges)}

function Base.IteratorEltype(::Type{T}) where {T <: SpeciesInteractionNetwork}
    return Base.HasEltype()
end

@testitem "SpeciesInteractionNetwork has both length and eltype" begin
    @test Base.IteratorSize(SpeciesInteractionNetwork{Bipartite,Binary}) == Base.HasLength()
    @test Base.IteratorEltype(SpeciesInteractionNetwork{Bipartite,Binary}) == Base.HasEltype()
end

function Base.isempty(N::T) where {T <: SpeciesInteractionNetwork}
    return iszero(length(N))
end

function _network_state(N::T, state::Integer) where {T <: SpeciesInteractionNetwork}
    C = findall(!iszero, N.edges.edges)[state]
    return (N.nodes[C]..., N[C])
end

function Base.iterate(N::T) where {T <: SpeciesInteractionNetwork}
    isempty(N) && return nothing
    return (_network_state(N, 1), 1)
end

function Base.iterate(N::T, state::Integer) where {T <: SpeciesInteractionNetwork}
    next = state == length(N) ? nothing : state+1
    isnothing(next) && return nothing
    return (_network_state(N, next), next)
end

@testitem "We can iterate over SpeciesInteractionNetwork" begin
    M = [true true true; true true false; false false true]
    edges = Binary(M)
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test [int for int in N][1] == (:A, :A, true)
    @test [int for int in N][length(N)] == (:C, :C, true)
    @test length([int for int in N]) == length(N)
end