"""
    distancetobase(::Type{SPM}, N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T, f) where {T, SPM <: ShortestPathMethod}

Measures the distance of species `sp` to a basal species in the food web, where
a basal species is defined as having a [`generaliry`](@ref) of 0.

Following *e.g.* [Thompson2012Food](@citet), we assign a distance to the base of
1 to the basal species. Primary consumers have a distance of 2, etc.

[Post2002long](@citet) notes that the trophic level can be obtained from the
maximum, mean, or minimum distance to a producer. Given that consumers may be
connected to more than one producer, one might argue that the mode or median of
these connections may be relevant. For this reason, the function `f` will
consume an array of distances, and return a scalar.

In favor of the minimum, one can argue that most energy transfer should happen
along short chains; but imagining a consumer atop a chain of length 5, also
connected directly to a producer, the minimum would give it a trophic level of
2, hiding its position at the top of the food web.

In favor of the maximum, one can argue that the higher chains give a better idea
of how far energy coming from the bottom of the food web can go. This is a
strong indication of how *vertically* diverse it is
[Duffy2007functional](@cite).

###### References

[Duffy2007functional](@citet*)

[Post2002long](@citet*)

[Thompson2012Food](@citet*)
"""
function distancetobase(::Type{SPM}, N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T, f) where {T, SPM <: ShortestPathMethod}
    paths = shortestpath(SPM, N, sp)
    basal_species = filter(s -> iszero(generality(N, s)), species(N))
    connected_basal_species = filter(s -> s in keys(paths), basal_species)
    if isempty(connected_basal_species)
        return 1.0
    end
    distances = [paths[s] for s in connected_basal_species]
    return f(distances) + 1.0
end

"""
    distancetobase(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T) where {T}

Default measure of [`distancetobase`](@ref) using the [`BellmanFord`](@ref)
shortest paths and the maximum distance.
"""
distancetobase(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T) where {T} = distancetobase(BellmanFord, N, sp, maximum)

"""
    distancetobase(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T, f) where {T}

Default measure of [`distancetobase`](@ref) using the [`BellmanFord`](@ref)
shortest paths and the distance returned by `f`.
"""
distancetobase(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T, f) where {T} = distancetobase(BellmanFord, N, sp, f)
