function Base.union(U::T, V::T) where {T <: Bipartite}
    return Bipartite(U.top ∪ V.top, U.bottom ∪ V.bottom)
end

@testitem "We can get the union of a bipartite partition" begin
    U = Bipartite([:A, :B, :C], [:a, :b, :c])
    V = Bipartite([:C, :D], [:a, :c, :d, :e])
    UV = U ∪ V
    for t in [:A, :B, :C, :D]
        @test t in species(UV, 1)
    end
    for b in [:a, :b, :c, :d, :e]
        @test b in species(UV, 2)
    end
end

function Base.union(U::T, V::T) where {T <: Unipartite}
    return Unipartite(U.margin ∪ V.margin)
end

@testitem "We can get the union of a unipartite partition" begin
    U = Unipartite([:A, :B, :C])
    V = Unipartite([:C, :D])
    UV = U ∪ V
    for s in [:A, :B, :C, :D]
        @test s in species(UV)
    end
end

function Base.union(U::T, V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
    nodes = U.nodes ∪ V.nodes
    edges = Binary(zeros(Bool, (richness(nodes,1), richness(nodes,2))))
    UV = SpeciesInteractionNetwork(nodes, edges)
    for u in interactions(U)
        UV[u[1],u[2]] = u[3]
    end
    for v in interactions(V)
        UV[v[1],v[2]] = v[3]
    end
    return UV
end

@testitem "We can get the union of a network" begin
    Nu = Unipartite([:A, :B, :C])
    Eu = Binary(rand(Bool, (richness(Nu,1), richness(Nu,2))))
    Nv = Unipartite([:C, :D])
    Ev = Binary(rand(Bool, (richness(Nv,1), richness(Nv,2))))
    U = SpeciesInteractionNetwork(Nu, Eu)
    V = SpeciesInteractionNetwork(Nv, Ev)
    UV = U ∪ V
    for s in species(U)
        @test s in species(UV)
    end
    for s in species(V)
        @test s in species(UV)
    end
    for i in U
        @test i in interactions(UV)
    end
    for i in V
        @test i in interactions(UV)
    end
end