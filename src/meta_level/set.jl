Base.union(U::T, V::T) where {T <: Bipartite} = Bipartite(U.top ∪ V.top, U.bottom ∪ V.bottom)
Base.union(U::T, V::T) where {T <: Unipartite} = Unipartite(U.margin ∪ V.margin)

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

@testitem "We can get the union of a unipartite partition" begin
    U = Unipartite([:A, :B, :C])
    V = Unipartite([:C, :D])
    UV = U ∪ V
    for s in [:A, :B, :C, :D]
        @test s in species(UV)
    end
end

Base.intersect(U::T, V::T) where {T <: Bipartite} = Bipartite(U.top ∩ V.top, U.bottom ∩ V.bottom)
Base.intersect(U::T, V::T) where {T <: Unipartite} = Unipartite(U.margin ∩ V.margin)

@testitem "We can get the intersect of a bipartite partition" begin
    U = Bipartite([:A, :B, :C], [:a, :b, :c])
    V = Bipartite([:C, :D], [:a, :c, :d, :e])
    UV = U ∩ V
    for t in [:C]
        @test t in species(UV, 1)
    end
    for b in [:a, :c]
        @test b in species(UV, 2)
    end
end

@testitem "We can get the intersect of a unipartite partition" begin
    U = Unipartite([:A, :B, :C])
    V = Unipartite([:C, :D])
    UV = U ∩ V
    for s in [:C]
        @test s in species(UV)
    end
end

function Base.union(U::T, V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
    nodes = U.nodes ∪ V.nodes
    edges = Binary(zeros(Bool, (richness(nodes,1), richness(nodes,2))))
    UV = SpeciesInteractionNetwork(nodes, edges)
    for u in interactions(U) ∪ interactions(V)
        UV[u[1],u[2]] = u[3]
    end
    return UV
end

@testitem "We can get the union of two networks" begin
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

function Base.intersect(U::T, V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
    nodes = U.nodes ∩ V.nodes
    edges = Binary(zeros(Bool, (richness(nodes,1), richness(nodes,2))))
    UV = SpeciesInteractionNetwork(nodes, edges)
    for u in interactions(U) ∩ interactions(V)
        UV[u[1],u[2]] = u[3]
    end
    return UV
end

@testitem "We can get the intersect of two networks" begin
    Nu = Unipartite([:A, :B, :C])
    Eu = Binary(rand(Bool, (richness(Nu,1), richness(Nu,2))))
    Nv = Unipartite([:C, :D])
    Ev = Binary(rand(Bool, (richness(Nv,1), richness(Nv,2))))
    U = SpeciesInteractionNetwork(Nu, Eu)
    V = SpeciesInteractionNetwork(Nv, Ev)
    UV = U ∩ V
    for s in species(U) ∩ species(V)
        @test s in species(UV)
    end
    for i in interactions(U)∩interactions(V)
        @test i in interactions(UV)
    end
end
