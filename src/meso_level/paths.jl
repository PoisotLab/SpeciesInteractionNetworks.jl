"""
    ShortestPathMethod

The first argument of [`shortestpath`](@ref) is a sub-type of
`ShortestPathMethod` which specifies the algorithm to use.

In ecological networks, the weight of interactions typically measure how *easy*
it is to move frome one node to the next. For this reason, we apply
transformations to various interaction types.

In binary networks, the weights and the distance are the same, because all
interactions have a value of 1.

For quantitative networks, we follow the approach of
[Newman2001Scientific](@citet), where the distance of an interaction is the
inverse of interaction strength. Nevertheless, it may be a good idea to correct
the interaction strength before applying the shortest path search, for example
by using [`normalize`](@ref). This normalization is *not* done by default, and
has to be explicit. Note that this normalisation is not changing the relative
ordering of paths, but is making the distances a little more comparable to the
binary case.

For probabilistic networks, the distance of an interaction is defined as ``1 +
(1 - p)``, so that an interaction happening with probability 1 has a distance of
one.

###### References

[Newman2001Scientific](@citet*)
"""
abstract type ShortestPathMethod end

_path_distance(::Type{Binary}, w) = w
_path_distance(::Type{Quantitative}, w) = one(w) / w
_path_distance(::Type{Probabilistic}, w) = 2one(w) - w


"""
    BellmanFord

The Bellman-Ford algorithmm [Shimbel1955Structure](@cite) measures the distance
from a source node to all other reachable nodes in the network.

Bellman-Ford is a weighted measure, in which the interactions are represented as
costs of moving from a node to another.

###### References

[Shimbel1955Structure](@citet*)
"""
abstract type BellmanFord <: ShortestPathMethod end

"""
    Dijkstra

###### References

[Dijkstra1959note](@citet*)
"""
abstract type Dijkstra <: ShortestPathMethod end

"""
    shortestpath(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, sp::T)

Defaults to [`BellmanFord`](@ref) for the shortest path algorithm. See also
[`Dijkstra`](@ref).

Note that in order to work with [`pathbetween`](@ref), any overload of
[`shortestpath`](@ref) *must* accept an `include_paths` keyword argument, that
when `true` returns a *second* dictionary giving the predecessors of each
reached node.
"""
shortestpath(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, sp::T; kwargs...) where {T} = shortestpath(BellmanFord, N, sp; kwargs...)

"""
    normalize(N::SpeciesInteractionNetwork{<:Partiteness, <:Quantitative})

Returns a quantitative network in which the interactions are normalized so that
the average of interactions is one. Note that this excludes self-interactions.
"""
function normalize(N::SpeciesInteractionNetwork{<:Partiteness, <:Quantitative})
    m = mean([i[3] for i in interactions(N) if i[1] != i[2]])
    return N ./ m
end

"""
    pathbetween(::Type{ShortestPathMethod}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, source::T, target::T) where {T}

Returns the path between `source` and `target`. The result is given as a vector
of interactions, *i.e.* it gives the subset of the output of `interactions(N)`
going from `source` to `target`.
"""
function pathbetween(::Type{SPM}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, source::T, target::T) where {SPM <: ShortestPathMethod, T}
    @assert source in species(N)
    @assert target in species(N)
    _, pred = shortestpath(SPM, N, source; include_paths=true)
    if !(target in keys(paths))
        return Vector{eltype(N)}()
    end

    path = eltype(N)[]

    reached = target
    while reached != source
        through = pred[reached]
        push!(path, (through, reached, N[through, reached]))
        pop!(pred, reached)
        reached = through
    end
    return reverse(path)

end

"""
    pathbetween(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, source::T, target::T) where {T}

Returns the path between `source` and `target`, using [`BellmanFord`](@ref) as
the default algorithm.
"""
function pathbetween(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, source::T, target::T) where {T}
    return pathbetween(BellmanFord, N, source, target)
end