function generality(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Binary}, sp::T) where {T}
    return sum(N[sp,:])
end

function generality(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Probabilistic}, sp::T) where {T}
    return sum(N[sp,:])
end

function generality(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Quantitative}, sp::T) where {T}
    return count(!iszero, N[sp,:])
end

function generality(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Binary}, sp::T) where {T}
    if sp in species(N,2)
        return 0
    end
    return sum(N[sp,:])
end

function generality(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Probabilistic}, sp::T) where {T}
    if sp in species(N,2)
        return 0.0
    end
    return sum(N[sp,:])
end

function generality(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Binary}, sp::T) where {T}
    if sp in species(N,2)
        return 0
    end
    return count(!iszero, N[sp,:])
end

function vulnerability(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Binary}, sp::T) where {T}
    if sp in species(N,1)
        return 0
    end
    return sum(N[:,sp])
end

function vulnerability(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Probabilistic}, sp::T) where {T}
    if sp in species(N,1)
        return 0.0
    end
    return sum(N[:,sp])
end

function vulnerability(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Binary}, sp::T) where {T}
    if sp in species(N,1)
        return 0
    end
    return count(!iszero, N[:,sp])
end

function vulnerability(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Binary}, sp::T) where {T}
    return sum(N[:,sp])
end

function vulnerability(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Probabilistic}, sp::T) where {T}
    return sum(N[:,sp])
end

function vulnerability(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Quantitative}, sp::T) where {T}
    return count(!iszero, N[:,sp])
end

function degree(N::SpeciesInteractionNetwork{<:Unipartite{T}, <:Interactions}, sp::T) where {T}
    d = generality(N, sp) + vulnerability(N, sp)
    correction = iszero(N[sp,sp]) ? zero(eltype(d)) : one(eltype(d))
    return d - correction
end

function degree(N::SpeciesInteractionNetwork{<:Bipartite{T}, <:Interactions}, sp::T) where {T}
    return generality(N, sp) + vulnerability(N, sp)
end

"""
    generality(N::SpeciesInteractionNetwork)

Returns the generality, *i.e.* (expected) number of interactions established by
each species. Note that you can specificy a second argument which is a species
from the network, giving the generality of this species alone.

###### References

[Schoener1989Food](@citet*)
"""
generality(N::SpeciesInteractionNetwork) = Dict([sp => generality(N, sp) for sp in species(N)])

"""
    vulnerability(N::SpeciesInteractionNetwork)

Returns the vulnerability, *i.e.* (expected) number of interactions received by
each species. Note that you can specificy a second argument which is a species
from the network, giving the vulnerability of this species alone.

###### References

[Schoener1989Food](@citet*)
"""
vulnerability(N::SpeciesInteractionNetwork) = Dict([sp => vulnerability(N, sp) for sp in species(N)])

"""
    degree(N::SpeciesInteractionNetwork)

Returns the degree, *i.e.* (expected) number of interactions involving each
species. Note that you can specificy a second argument which is a species from
the network, giving the degree of this species alone.

###### References

[Schoener1989Food](@citet*)
"""
degree(N::SpeciesInteractionNetwork) = Dict([sp => degree(N, sp) for sp in species(N)])

@testitem "We can get the degree of a unipartite network with self-loops" begin
    nodes = Unipartite([:A, :B, :C, :D])
    edges = Binary(Bool[1 1 1 1; 0 0 0 1; 0 0 1 1; 1 0 0 1])
    N = SpeciesInteractionNetwork(nodes, edges)
    D = degree(N)
    @test D[:A] == 5
    @test D[:B] == 2
    @test D[:C] == 3
    @test D[:D] == 5
end

@testitem "We can get the degree using the species-specific call" begin
    nodes = Unipartite([:A, :B, :C, :D])
    edges = Binary(Bool[1 1 1 1; 0 0 0 1; 0 0 1 1; 1 0 0 1])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test degree(N, :A) == 5
end