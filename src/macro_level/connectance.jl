"""
    links(N::SpeciesInteractionNetwork{<:Partiteness, <:Interactions})

The number of links in a network is defined as the number of non-zero elements
in the network. For quantitative networks, this is *not* weighted by the
intensity of each interaction.
"""
function links(N::SpeciesInteractionNetwork)
    return count(!iszero, N.edges.edges)
end

"""
    links(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})

The number of links in a probabilistic network is defined as the expected number
of interactions in the network. See also [`links_variance`](@ref).

###### References

[Poisot2015structure](@citet*)
"""
function links(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})
    return sum(N.edges.edges)
end

"""
    connectance(N::SpeciesInteractionNetwork)

The connectance, also known as the density, of a network is defined as the ratio
between the number of existing links and the number of possible links. This
value is in the unit interval. For a probabilistic network, this returns the
*expected* connectance.

###### References

[Martinez1992Constant](@citet*)
"""
function connectance(N::SpeciesInteractionNetwork)
    return links(N) / (richness(N,1)*richness(N,2))
end

"""
    connectance(N::SpeciesInteractionNetwork)

The linkage density of a network is the average number of interaction per
species (or its expected value in the case of probabilistic networks).
"""
function linkagedensity(N::SpeciesInteractionNetwork)
    return links(N) / richness(N)
end

"""
    links_variance(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})

The variance in the expected number of links in a probabilistic network is
defined as the sum of ``p \\times (1 - p)``, where ``p`` is the probability of
all interactions (including ``p = 0``).

###### References

[Poisot2015structure](@citet*)
"""
function links_variance(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})
    return links(N .* (1 .- N))
end


"""
    connectance_variance(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})

See also [`connectance`](@ref).

###### References

[Poisot2015structure](@citet*)
"""
function connectance_variance(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})
    return links_variance(N) / (richness(N,1)*richness(N,2))
end

"""
    linkagedensity_variance(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})

See also [`linkagedensity_variance`](@ref).

###### References

[Poisot2015structure](@citet*)
"""
function linkagedensity_variance(N::SpeciesInteractionNetwork{<:Partiteness, <:Probabilistic})
    return links_variance(N) / (richness(N))
end