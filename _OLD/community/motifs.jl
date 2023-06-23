function _inner_find_motif(N::T1, m::T2) where {T1<:UnipartiteNetwork,T2<:UnipartiteNetwork}
    motif_permutations = _permute_motif(m)
    matching_species = []
    for species_combination in combinations(1:richness(N), richness(m))
        if N.edges[species_combination,species_combination] ∈ motif_permutations
            push!(matching_species, (species_combination,))
        end
    end
    return matching_species
end

function _inner_find_motif(N::T1, m::T2) where {T1<:BipartiteNetwork,T2<:BipartiteNetwork}
    motif_permutations = _permute_motif(m)
    matching_species = []
    top_combinations = combinations(species(N; dims=1), richness(m; dims=1))
    bottom_combinations = combinations(species(N; dims=2), richness(m; dims=2))
    # Pre-allocate the species positions
    top_sp_pos = zeros(Int64, richness(m; dims=1))
    bot_sp_pos = zeros(Int64, richness(m; dims=2))
    # Positions of species (they don't change!)
    p1 = Dict(zip(species(N; dims=1), 1:richness(N; dims=1)))
    p2 = Dict(zip(species(N; dims=2), 1:richness(N; dims=2)))
    for top_species in top_combinations
        for i in 1:length(top_species)
            top_sp_pos[i] = p1[top_species[i]]
        end
        for bottom_species in bottom_combinations
            for i in 1:length(bottom_species)
                bot_sp_pos[i] = p2[bottom_species[i]]
            end
            if N.edges[top_sp_pos, bot_sp_pos] ∈ motif_permutations
                push!(matching_species, (top_species, bottom_species))
            end
        end
    end
    return matching_species
end

function _inner_find_motif(N::T1, m::T2) where {T1<:UnipartiteProbabilisticNetwork, T2<:UnipartiteNetwork}
    motif_permutations = _permute_motif(m)
    all_combinations = []
    for species_combination in combinations(species(N), richness(m))
        isg = N[species_combination]
        motif_mean, motif_var = 0.0, 0.0
        for perm in motif_permutations
            imat = zeros(eltype(N.edges), size(m))
            for i in eachindex(imat)
                imat[i] = (perm[i] ? 2.0*isg[i] : 1.0)-isg[i]
            end
            pmm, pmv = prod(imat), EcologicalNetworks._multiplicative_bernoulli_variance(imat)
            motif_mean += pmm
            motif_var += pmv
        end
        push!(all_combinations, ((species_combination,),(motif_mean, motif_var)))
    end
    return all_combinations
end

function _inner_find_motif(N::T1, m::T2) where {T1<:BipartiteProbabilisticNetwork, T2<:BipartiteNetwork}
    motif_permutations = _permute_motif(m)
    all_combinations = []
    top_combinations = combinations(species(N; dims=1), richness(m; dims=1))
    bottom_combinations = combinations(species(N; dims=2), richness(m; dims=2))
    for top_species in top_combinations, bottom_species in bottom_combinations
        isg = N[top_species, bottom_species]
        motif_mean, motif_var = 0.0, 0.0
        for perm in motif_permutations
            imat = zeros(eltype(N.edges), size(m))
            for i in eachindex(imat)
                imat[i] = (perm[i] ? 2.0*isg[i] : 1.0)-isg[i]
            end
            pmm, pmv = prod(imat), EcologicalNetworks._multiplicative_bernoulli_variance(imat)
            motif_mean += pmm
            motif_var += pmv
        end
        push!(all_combinations, ((top_species, bottom_species),(motif_mean, motif_var)))
    end
    return all_combinations
end

"""
    find_motif(N::T1, m::T2) where {T1<:AbstractEcologicalNetwork, T2<:BinaryNetwork}

Returns an array of tuples, in which each tuple contains the species that are
part of the motif. The length of the array gives the number of times the motif
was found. For probabilistic networks, the tuple also contains the probability
of observing the species in the correct conformation for the motif, as well as
the variance.

#### References

- Milo, R., Shen-Orr, S., Itzkovitz, S., Kashtan, N., Chklovskii, D., Alon, U.,
  2002. Network motifs: simple building blocks of complex networks. Science 298,
  824–7. https://doi.org/10.1126/science.298.5594.824

- Poisot, T., Cirtwill, A.R., Cazelles, K., Gravel, D., Fortin, M.-J., Stouffer,
  D.B., 2016. The structure of probabilistic networks. Methods in Ecology and
  Evolution 7, 303–312. https://doi.org/10.1111/2041-210X.12468
"""
function find_motif(N::T1, m::T2) where {T1<:AbstractEcologicalNetwork, T2<:BinaryNetwork}
    M = copy(N)
    if typeof(N) <: AbstractUnipartiteNetwork
        M = nodiagonal(M)
    end
    if typeof(M) <: QuantitativeNetwork
        M = convert(BinaryNetwork, M)
    end
    return _inner_find_motif(M, m)
end

"""
    expected_motif_count(s)

Get the expected number of motifs (and variance) from the output of `find_motif`
on a probabilistic network.

#### References

- Poisot, T., Cirtwill, A.R., Cazelles, K., Gravel, D., Fortin, M.-J., Stouffer,
  D.B., 2016. The structure of probabilistic networks. Methods in Ecology and
  Evolution 7, 303–312. https://doi.org/10.1111/2041-210X.12468
"""
function expected_motif_count(s)
    m = [x[2][1] for x in s]
    v = [x[2][2] for x in s]
    return (sum(m), sum(v))
end
