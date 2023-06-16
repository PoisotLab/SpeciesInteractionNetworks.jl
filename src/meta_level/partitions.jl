"""
    BetaDivComponent

The [`betadiversity`](@ref) methods all use a subtype of `BetaDivComponent` as
their first argument, to determine which component should be measured.

All of the partitions follow the [Koleff2003Measuring](@cite) approach, where
the beta diversity is measured on the cardinality of sets.
"""
abstract type BetaDivComponent end

"""
    βS

Specie
"""
abstract type βS <: BetaDivComponent end

"""
    βOS
"""
abstract type βOS <: BetaDivComponent end

"""
    βWN

"""
abstract type βWN <: BetaDivComponent end

"""
    betadiversity(::Type{βS},U::T,V::T,dims::Integer = 0,) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}

Species-level β diversity between networks. By default, this return the species
dissimilarity for the entire network. An optional last argument `dims` can be
used, to specify top (`1`) and bottom (`2`) levels.

###### References

[Koleff2003Measuring](@cite)
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

"""
    betadiversity(::Type{βWN},U::T,V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}

Network-level β diversity between networks. By default, this return the species
dissimilarity for the entire network. An optional last argument `dims` can be
used, to specify top (`1`) and bottom (`2`) levels.

###### References

[Poisot2012dissimilarity](@cite)

[Canard2014Empirical](@cite)

[Poisot2022Dissimilarity](@cite)
"""
function betadiversity(
    ::Type{βWN},
    U::T,
    V::T
) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
    shared = length(intersect(U, V))
    right = length(V) - shared
    left = length(U) - shared
    return (; left, shared, right)
end

@testitem "We can get the βWN for a unipartite network" begin
    Nu = Unipartite([:A, :B, :C])
    Eu = Binary(zeros(Bool, (richness(Nu,1), richness(Nu,2))))
    Nv = Unipartite([:A, :C, :D])
    Ev = Binary(zeros(Bool, (richness(Nv,1), richness(Nv,2))))
    U = SpeciesInteractionNetwork(Nu, Eu)
    V = SpeciesInteractionNetwork(Nv, Ev)
    U[:A, :C] = true
    U[:A, :B] = true
    U[:B, :C] = true
    V[:A, :C] = true
    V[:C, :D] = true
    V[:A, :D] = true
    part = betadiversity(βWN, U, V)
    @test part.shared == 1
    @test part.right == 2
    @test part.left == 2
end

function betadiversity(
    ::Type{βOS},
    U::T,
    V::T
) where {T <: SpeciesInteractionNetwork{<:Unipartite, <:Binary}}
    core = species(U∩V)
    if isempty(core)
        left = length(U)
        right = length(V)
        shared = 0
        return (; left, shared, right)
    end
    Us = subgraph(U, core)
    Vs = subgraph(V, core)
    return betadiversity(βWN, Us, Vs)
end

function betadiversity(
    ::Type{βOS},
    U::T,
    V::T
) where {T <: SpeciesInteractionNetwork{<:Bipartite, <:Binary}}
    core_top = species(U∩V, 1)
    core_bottom = species(U∩V, 2)
    if isempty(core_top)|isempty(core_bottom)
        left = length(U)
        right = length(V)
        shared = 0
        return (; left, shared, right)
    end
    Us = subgraph(U, core_top, core_bottom)
    Vs = subgraph(V, core_top, core_bottom)
    return betadiversity(βWN, Us, Vs)
end

@testitem "We can get the βOS for a unipartite network" begin
    Nu = Unipartite([:A, :B, :C])
    Eu = Binary(zeros(Bool, (richness(Nu,1), richness(Nu,2))))
    Nv = Unipartite([:A, :C, :D])
    Ev = Binary(zeros(Bool, (richness(Nv,1), richness(Nv,2))))
    U = SpeciesInteractionNetwork(Nu, Eu)
    V = SpeciesInteractionNetwork(Nv, Ev)
    U[:A, :C] = true
    U[:A, :B] = true
    U[:B, :C] = true
    V[:A, :C] = true
    V[:C, :D] = true
    V[:A, :D] = true
    part = betadiversity(βOS, U, V)
    @test part.shared = 1
    @test part.left = 0
    @test part.right = 0
end

@testitem "We get the correct βOS when no species are shared" begin
    Nu = Unipartite([:A, :B, :C])
    Eu = Binary(zeros(Bool, (richness(Nu,1), richness(Nu,2))))
    Nv = Unipartite([:D, :E, :F])
    Ev = Binary(zeros(Bool, (richness(Nv,1), richness(Nv,2))))
    U = SpeciesInteractionNetwork(Nu, Eu)
    V = SpeciesInteractionNetwork(Nv, Ev)
    part = betadiversity(βOS, U, V)
    @test part.shared = 0
    @test part.left = length(U)
    @test part.right = length(V)
end