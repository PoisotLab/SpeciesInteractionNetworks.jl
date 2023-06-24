"""
    trophic_level(N::UnipartiteNetwork; kwargs...)

Returns the *fractional* trophic level (after Odum & Heald 1975)
of species in a binary unipartite network. The trophic level is
calculated after Pauly & Christensen (1995), specifically as TLᵢ = 1 +
∑ⱼ(TLⱼ×DCᵢⱼ)/∑ⱼ(DCᵢⱼ), wherein TLᵢ is the trophic level
of species i, and DCᵢⱼ is the proportion of species j in the diet of
species i. Note that this function is calculated on a network where the
diagonal (i.e. self loops) are removed.

This function uses a *pseudo*-inverse to solve the linear system
described above *iff* the determinant of the diet matrix is 0 (it is
non-invertible). Otherwise, it will use an exact inverse.
"""
function trophic_level(N::UnipartiteNetwork)
    # Get the degree as a vector, ordered the same way as species
    kout = EcologicalNetworks.degree_out(N)
    𝐤 = [kout[s] for s in species(N)]
    # Adjacency matrix to solve the TL
    𝐀 = adjacency(N)
    # Diet matrix
    𝐃 = .-(𝐀 ./ 𝐤)
    replace!(𝐃, NaN => -0.0)
    𝐃[diagind(𝐃)] .= 1.0 .- 𝐃[diagind(𝐃)]
    # Solve with the inverse matrix
    𝐛 = ones(eltype(𝐃), richness(N))
    invfunc = iszero(det(𝐃)) ? pinv : inv
    return Dict(zip(species(N), invfunc(𝐃) * 𝐛))
end

"""
    omnivory(N::T) where {T <: UnipartiteNetwork}

Returns a vector of length richness(N) with an index of consumption across
trophic levels. Note we explicitly assume that diet contribution is spread
evenly across prey species, in proportion to their degree. Omnivory is
measured based on the output of `trophic_level`.

#### References
- Christensen, Villy, and Daniel Pauly. "ECOPATH II—a software for balancing
  steady-state ecosystem models and calculating network characteristics."
  Ecological modelling 61, no. 3-4 (1992): 169-185.

"""
function omnivory(N::T) where {T<:UnipartiteNetwork}
    # Get the trophic level as an array
    TL = trophic_level(N)

    # Prepare the omnivory index dictionary
    OI = Dict(zip(species(N), fill(0.0, richness(N))))

    # Loop
    for sp in species(N)
        if !isempty(N[sp,:])
            TLj = [TL[prey] for prey in N[sp,:]]
            # Note that we don't need the normalisation here, because ∑Dᵢⱼ=1
            OI[sp] = sum((TLj .- mean(TLj)).^2.0)
        end
    end

    return OI
end
