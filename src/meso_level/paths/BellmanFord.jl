_bf_distance(::Type{Binary}, w) = w
_bf_distance(::Type{Quantitative}, w) = one(w) / w
_bf_distance(::Type{Probabilistic}, w) = 2one(w) - w

"""
    normalize(N::SpeciesInteractionNetwork{<:Partiteness, <:Quantitative})

Returns a quantitative network in which the interactions are normalized so that
the average of interactions is one. Note that this excludes self-interactions.
"""
function normalize(N::SpeciesInteractionNetwork{<:Partiteness, <:Quantitative})
    m = mean([i[3] for i in interactions(N) if i[1] != i[2]])
    return N ./ m
end

function shortestpath(
    ::Type{BellmanFord},
    N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions},
    sp::T,
) where {T}
    @assert sp in species(N)

    dist = Dict([s => Inf for s in species(N)])
    pred = Dict{T, Union{Nothing, T}}([s => nothing for s in species(N)])
    dist[sp] = 0.0

    pool = interactions(N)

    changes_made = true

    df = (x) -> _bf_distance(_edgetype(N), x)

    for _ in length(N) - 1
        changes_made = false
        for int in pool
            if (dist[int[1]] + df(int[3])) < dist[int[2]]
                dist[int[2]] = dist[int[1]] + df(int[3])
                pred[int[2]] = int[1]
                changes_made = true
            end
        end
        if !changes_made
            break
        end
    end

    for s in species(N)
        if (isinf(dist[s])) | (isnothing(pred[s]))
            pop!(dist, s, nothing)
        end
    end
    return dist
end

@testitem "Bellman-Ford works on binary networks" begin
    nodes = Unipartite([:A, :B, :C, :D, :E, :F, :G])
    edges = Binary(zeros(Bool, (richness(nodes), richness(nodes))))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in
        [(:A, :B), (:B, :C), (:C, :D), (:B, :E), (:C, :F), (:E, :F), (:F, :G), (:D, :D)]
        N[edge...] = true
    end
    bf = shortestpath(BellmanFord, N, :B)
    @test bf[:C] == 1
    @test bf[:D] == 2
    @test bf[:E] == 1
    @test bf[:F] == 2
    @test bf[:G] == 3
end

@testitem "Bellman-Ford works on quantitative networks" begin
    nodes = Unipartite([:A, :B, :C, :D, :E, :F, :G])
    edges = Quantitative(zeros(Int8, (richness(nodes), richness(nodes))))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in [
        (:A, :B, 1),
        (:B, :C, 1),
        (:C, :D, 2),
        (:B, :E, 2),
        (:C, :F, 1),
        (:E, :F, 3),
        (:F, :G, 4),
        (:D, :D, 9),
        (:B, :A, 4),
    ]
        N[edge[1], edge[2]] = edge[3]
    end
    bf = shortestpath(BellmanFord, normalize(N), :B)
    @test bf[:A] == 0.5625
    @test bf[:C] == 2.25
    @test bf[:D] == 3.375
    @test bf[:E] == 1.125
    @test bf[:F] == 1.875
    @test bf[:G] == 2.4375
end

@testitem "Bellman-Ford works on probabilistic networks" begin
    nodes = Unipartite([:A, :B, :C, :D])
    edges = Probabilistic(zeros(Float16, (richness(nodes), richness(nodes))))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in
        [(:A, :B, 0.8), (:A, :C, 0.5), (:A, :D, 0.1), (:B, :D, 1 / 3), (:C, :D, 0.9)]
        N[edge[1], edge[2]] = edge[3]
    end
    bf = shortestpath(BellmanFord, N, :A)
    @test bf[:C] < bf[:D]
    @test bf[:B] < bf[:D]
    @test bf[:B] < bf[:C]
end