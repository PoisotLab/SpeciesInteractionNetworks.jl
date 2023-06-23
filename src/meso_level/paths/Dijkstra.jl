function shortestpath(::Type{Dijkstra}, N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, sp::T) where {T}
    @assert sp in species(N)

    dist = Dict([s => Inf for s in species(N)])
    pred = Dict{T, Union{Nothing, T}}([s => nothing for s in species(N)])
    dist[sp] = 0.0
    Q = species(N)

    df = (x) -> SpeciesInteractionNetworks._path_distance(SpeciesInteractionNetworks._edgetype(N), x)

    while !isempty(Q)
        _, u = findmin(filter(p -> p.first ∈ Q, dist))
        setdiff!(Q, [u])
        whereto = filter(v -> v ∈ Q, successors(N, u))
        if isempty(whereto)
            break
        end
        for v in whereto
            proposal = dist[u] + df(N[u,v])
            if proposal < dist[v]
                dist[v] = proposal
                pred[v] = u
            end
        end
    end

    for s in species(N)
        if (isinf(dist[s])) | (isnothing(pred[s]))
            pop!(dist, s, nothing)
        end
    end
    return dist
end

@testitem "Dijkstra works on binary networks" begin
    nodes = Unipartite([:A, :B, :C, :D, :E, :F, :G])
    edges = Binary(zeros(Bool, (richness(nodes), richness(nodes))))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in
        [(:A, :B), (:B, :C), (:C, :D), (:B, :E), (:C, :F), (:E, :F), (:F, :G), (:D, :D)]
        N[edge...] = true
    end
    bf = shortestpath(Dijkstra, N, :B)
    @test bf[:C] == 1
    @test bf[:D] == 2
    @test bf[:E] == 1
    @test bf[:F] == 2
    @test bf[:G] == 3
end