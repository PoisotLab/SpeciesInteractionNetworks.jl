abstract type ShortestPathMethod end

"""
    BellmanFord

The Bellman-Ford algorithmm [Shimbel1955Structure](@cite) measures the distance
from a source node to all other reachable nodes in the network.

Bellman-Ford is a weighted measure, in which the interactions are represented as
costs of moving from a node to another. In ecological networks, the weight of
interactions typically measure how *easy* it is to move frome one node to the
next. For this reason, we apply transformations to various interaction types.

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
1.

###### References

[Newman2001Scientific](@citet*)

[Shimbel1955Structure](@citet*)
"""
abstract type BellmanFord <: ShortestPathMethod end

abstract type Djikstra <: ShortestPathMethod end

shortestpath(N::SpeciesInteractionNetwork) = path(BellmanFord, N)

#=
function bellman_ford(N::T, source::K) where {T <: DeterministicNetwork, K}

    source ∈ species(N) || throw(ArgumentError("Species $(source) is not part of the network"))

    d = Dict([s => Inf64 for s in species(N)])

    # The dictionary for the predecessor start as empty, and this saves some
    # issues with the species types being multiple possible types
    p = Dict{K,K}()
    # We will sizehint! it for good measure, but it may not be entirely filled
    sizehint!(p, richness(N))

    d[source] = 0.0

    all_edges = interactions(N)

    for i in 1:(richness(N)-1)
        for interaction in all_edges
            w = get(interaction, :strength, 1.0)
            if d[interaction.to] > (d[interaction.from] + w)
                d[interaction.to] = d[interaction.from] + w
                p[interaction.to] = interaction.from
            end
        end
    end

    return [(from=source, to=k, weight=d[k]) for (k,v) in p]
end


"""
    bellman_ford(N::T) where {T <: DeterministicNetwork}

Bellman-ford algorithm to return the shortest / easiest paths between all pairs
of species in the networks, as long as paths exists. This function will return a
tuple, with fields `from`, `to` and `weight`. The number of elements in the
tuple is the number of paths. This function works with quantitative and binary
networks, and assumes that no interactions are negative.

Currently, the Bellman-Ford algorithm is *slower* than the `shortest_path`
function, but the arguments are returned in a more usable way. Note that the
speed penalty is only valid when measuring the shortest paths in the entire
network (and will be fixed relatively soon), and does not apply as much for the
shortest paths from a single source node.
"""
function bellman_ford(N::T) where {T <: DeterministicNetwork}
    global paths
    @inbounds for i in 1:richness(N)
        i == 1 && (paths = bellman_ford(N, species(N)[i]))
        i == 1 || append!(paths, bellman_ford(N, species(N)[i]))
    end
    return paths
end

function get_adj_list(N::T, species::Array{K,1}) where {T <: DeterministicNetwork, K}
    adj_list = Dict{K,Array{Tuple{Float64,K}}}()
    for s in species
        adj_list[s] = []
    end
    for interaction in interactions(N)
        w = get(interaction, :strength, 1.0)
        push!(adj_list[interaction.from], (w, interaction.to))
    end
    return adj_list
end

"""
    dijkstra(N::T) where {T <: DeterministicNetwork}

Dijkstra algorithm to return the shortest / easiest paths between all pairs
of species in the networks, as long as paths exists. This function will return a
tuple, with fields `from`, `to` and `weight`. The number of elements in the
tuple is the number of paths. This function works with quantitative and binary
networks, and assumes that no interactions are negative.
"""
function dijkstra(N::T) where {T <: DeterministicNetwork}
    #QUESTION dealing with bipartite networks?

    species_of_N = species(N)
    K = eltype(species_of_N)

    d = Dict([(s1, s2) => Inf64 for s1 in species_of_N for s2 in species_of_N])

    p = Dict{Tuple{K,K},K}()

    adj_list = get_adj_list(N, species_of_N)

    to_check = Set{K}()
    sizehint!(to_check, richness(N))

    for s in species_of_N
        d[(s,s)] = 0.0
        push!(to_check, s)
    end

    while length(to_check) > 0
        current = pop!(to_check)
        for (w_cur2neig, neighbor) in adj_list[current]
            # check if there is a shorter path from node to neighbor via current
            for s in species_of_N
                if d[s, neighbor] > (d[s, current] + w_cur2neig)
                    d[s, neighbor] = d[s, current] + w_cur2neig
                    p[(s, neighbor)] = current
                    push!(to_check, neighbor)
                end
            end
        end
    end
    return [(from=s, to=n, weight=d[(s,n)]) for ((s,n),c) in p]
end


"""
    dijkstra(N::T, source::K) where {T <: DeterministicNetwork, K}

Dijkstra's algorithm to return the shortest / easiest paths from a source
species. Refer to the `bellman_ford` global documentation for the output format.
"""
function dijkstra(N::T, source::K) where {T <: DeterministicNetwork, K}

    source ∈ species(N) || throw(ArgumentError("Species $(source) is not part of the network"))

    d = Dict([s => Inf64 for s in species(N)])

    # The dictionary for the predecessor start as empty, and this saves some
    # issues with the species types being multiple possible types
    p = Dict{K,K}()
    # We will sizehint! it for good measure, but it may not be entirely filled
    sizehint!(p, richness(N))

    d[source] = 0.0

    adj_list = get_adj_list(N, species(N))

    #to_check = PriorityQueue{K,Float64}()
    #append!(to_check, (0.0, source))

    to_check = [(0.0, source)]

    while length(to_check) > 0
        dist_source_current, current = heappop!(to_check)
        for (w, neighbor) in adj_list[current]  # scan neighbors
            dist_via_neighbor = (dist_source_current + w)
            if d[neighbor] > dist_via_neighbor
                d[neighbor] = dist_via_neighbor
                p[neighbor] = current
                if haskey(adj_list, neighbor)
                    heappush!(to_check, (dist_via_neighbor, neighbor))
                end
            end
        end
    end

    return [(from=source, to=k, weight=d[k]) for (k,v) in p]
end
=#