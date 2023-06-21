abstract type ShortestPathMethod end
abstract type BellmanFord <: ShortestPathMethod end
abstract type Djikstra <: ShortestPathMethod end

shortestpath(N::SpeciesInteractionNetwork) = path(BellmanFord, N)

_bf_distance(d::T, ::Binary{T}) where {T} = d
_bf_distance(d::T, ::Probabilistic{T}) where {T} = 2.0 - d
_bf_distance(d::T, ::Quantitative{T}) where {T} = d ^(-2.0)

"""
    shortestpath(::Type{BellmanFord}, N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary,Quantitative}}, source)

Uses the Bellman-Ford algorithm ...

For binary networks, all interactions incur a distance of one. For quantitative
networks, an interaction with edge weight w incurs a distance of w⁻², so that
strong interactions pull nodes together.
"""
function shortestpath(::Type{BellmanFord}, N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary,Quantitative}}, source)
    
    distance = fill(Inf, richness(N))
    pred = Vector{eltype(N.nodes)}(undef, richness(N))
    
    source_id = findfirst(isequal(source), species(N))
    distance[source_id] = 0.0
    
    for _ in 1:(richness(N)-1)
        changes_made = false
        for interaction in N
            from = findfirst(isequal(interaction[1]), species(N))
            to = findfirst(isequal(interaction[2]), species(N))
            if (distance[from] + interaction[3]) < distance[to]
                distance[to] = distance[from] + _bf_distance(interaction[3], N.edges)
                pred[to] = interaction[1]
                changes_made = true
            end
        end
        if !(changes_made)
            break
        end
    end
    
    found_paths = unique([(species(N)[p],distance[p]) for p in findall(d -> !(iszero(d)|isinf(d)), distance)])
    return Dict(zip(first.(found_paths), last.(found_paths)))

end

@testitem "Bellman-Ford works on binary networks" begin
    nodes = Unipartite([:a, :b, :c, :d, :e])
    edges = Binary(zeros(Bool, (richness(nodes), richness(nodes))))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in [(:a,:b), (:b,:c), (:c,:d), (:c,:e)]
        N[edge...] = true
    end
    bf = shortestpath(BellmanFord, N, :a)
    @test bf[:b] == 1.0
    @test bf[:c] == 2.0
    @test bf[:d] == bf[:e] == 3.0
    @test !(:a in keys(bf))
end

@testitem "Bellman-Ford works on quantitative networks" begin
    nodes = Unipartite([:a, :b, :c, :d, :e, :f])
    edges = Quantitative(fill(0, (richness(nodes), richness(nodes))))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in [(:a,:b,1), (:b,:c,1), (:c,:d,1), (:c,:e,2), (:e,:f,1), (:d,:f,1)]
        N[edge[1:2]...] = edge[3]
    end
    bf = shortestpath(BellmanFord, N, :a)
    @test bf[:b] == 1.0
    @test bf[:c] == 2.0
    @test bf[:d] == 3.0
    @test bf[:e] < 3.0
    @test bf[:e] == 2 + 2^-2
    @test bf[:f] == 3 + 2^-2
    @test !(:a in keys(bf))
end

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