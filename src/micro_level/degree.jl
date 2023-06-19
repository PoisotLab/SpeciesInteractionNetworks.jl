"""
    degree(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

The degree of a species in a unipartite network is the sum of its generality and
vulnerability. This function *will* count self-links as part of the degree.
See also [`generality`](@ref), [`vulnerability`](@ref).
"""
function degree(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    deg = Dict([sp => length(predecessors(N, sp))+length(successors(N, sp)) for sp in species(N)])
    for sp in species(N)
        if N[sp,sp]
            deg[sp] -= 1
        end
    end
    return deg
end

"""
    vulnerability(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

The in-degree of a species in a unipartite network is the number of its
successors. This function *will* count self-links as part of the vulnerability.
See also [`degree`](@ref), [`generality`](@ref).

Note that by contrast to the original definition of vulnerability
[Schoener1989Food](@cite), this version is *not* corrected for network richness.

###### References

[Schoener1989Food](@citet)
"""
function vulnerability(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    deg = Dict([sp => length(predecessors(N, sp)) for sp in species(N)])
    return deg
end

"""
    generality(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

The in-degree of a species in a unipartite network is the number of its
successors. This function *will* count self-links as part of the generality. See
also [`degree`](@ref), [`vulnerability`](@ref).

Note that by contrast to the original definition of vulnerability
[Schoener1989Food](@cite), this version is *not* corrected for network richness.

###### References

[Schoener1989Food](@citet)
"""
function generality(N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    deg = Dict([sp => length(successors(N, sp)) for sp in species(N)])
    return deg
end

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

"""
    vulnerability(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})

In a bipartite network, the vulnerability is only defined for bottom-level
species, and is the number of their predecessors (as given by the
[`predecessors`](@ref) function).
"""
function vulnerability(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})
    return Dict([sp => length(predecessors(N, sp)) for sp in species(N,2)])
end

"""
    generality(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})

In a bipartite network, the generality is only defined for top-level species,
and is the number of their successors (as given by the [`successors`](@ref)
function).
"""
function generality(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})
    return Dict([sp => length(successors(N, sp)) for sp in species(N,1)])
end

"""
    degree(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})

In a bipartite network, the degree is the combination of the
[`generality`](@ref) and [`vulnerability`](@ref) of all species.
"""
function degree(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})
    return merge(generality(N), vulnerability(N))
end