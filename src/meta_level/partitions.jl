"""
    BetaDivComponent

The [`betadiversity`](@ref) methods all use a subtype of `BetaDivComponent` as
their first argument, to determine which component should be measured.
"""
abstract type BetaDivComponent end

"""
    βS

Specie
"""
abstract type βS <: BetaDivComponent end

"""
    βOS

Specie
"""
abstract type βOS <: BetaDivComponent end

"""
    βWN

Specie
"""
abstract type βWN <: BetaDivComponent end

"""
    betadiversity(::Type{βS},U::T,V::T,dims::Integer = 0,) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}

Species-level β diversity between networks. By default, this return the species
dissimilarity for the entire network. An optional last argument `dims` can be
used, to specify top (`1`) and bottom (`2`) levels.
"""
function betadiversity(
    ::Type{βS},
    U::T,
    V::T,
    dims::Integer = 0,
) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
    shared = richness(intersect(U, V), dims)
    right = richness(V, dims) - shared
    left = richness(U, dims) - shared
    return (; left, shared, right)
end

@testitem "We can get the βS for a unipartite network" begin
    Nu = Unipartite([:A, :B, :C])
    Eu = Binary(rand(Bool, (richness(Nu,1), richness(Nu,2))))
    Nv = Unipartite([:C, :D])
    Ev = Binary(rand(Bool, (richness(Nv,1), richness(Nv,2))))
    U = SpeciesInteractionNetwork(Nu, Eu)
    V = SpeciesInteractionNetwork(Nv, Ev)
    part = betadiversity(βS, U, V)
    @test part.shared == 1
    @test part.right == 1
    @test part.left == 2
end 