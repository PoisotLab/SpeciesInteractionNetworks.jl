Tables.istable(::Type{<:SpeciesInteractionNetwork}) = true

Tables.rowaccess(::Type{<:SpeciesInteractionNetwork}) = true
Tables.rows(N::SpeciesInteractionNetwork) = interactions(N)

Tables.columnaccess(::Type{<:SpeciesInteractionNetwork}) = true
Tables.getcolumn(N::SpeciesInteractionNetwork, i::Int) = [interaction[i] for interaction in interactions(N)]
Tables.getcolumn(N::SpeciesInteractionNetwork, nm::Symbol) = [interaction[findfirst(isequal(nm), Tables.columnnames(N))] for interaction in interactions(N)]
Tables.columnnames(::SpeciesInteractionNetwork) = (:from, :to, :value)

Tables.schema(::SpeciesInteractionNetwork{<:Partiteness{T1}, <:Interactions{T2}}) where {T1, T2} = Tables.Schema([:from, :to, :value], [T1, T1, T2])

@testitem "The networks behave as tables" begin
    import SpeciesInteractionNetworks.Tables
    M = [true true true; true true false; false false true]
    edges = Binary(M)
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test Tables.istable(typeof(N))
    @test Tables.rowaccess(typeof(N))
    @test Tables.rows(N) == interactions(N)
    @test Tables.columnaccess(typeof(N))
    @test Tables.columnnames(N) == (:from, :to, :value)
end