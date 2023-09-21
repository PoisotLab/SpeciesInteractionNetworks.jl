# # The niche model of food webs

#=
!!! example "The niche model of food webs"
    
    In this example, we will look at data from 50 pelagic food webs
    [Havens1992Scale](@cite), to figure out how the niche model of food webs
    captures the structure of empirical networks.
=#

using SpeciesInteractionNetworks
import Statistics
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide

# In order to carry out this example, we will explicitely import functions from
# the `Mangal` package, which will allow to query the
# [`mangal.io`](https://mangal.io/) database. This is done with
# [`mangalnetwork`](@ref).

import SpeciesInteractionNetworks.Mangal
dataset = Mangal.dataset("havens_1992")
networks = mangalnetwork.(Mangal.networks(dataset); taxonlevel=true)

@info extrema(richness.(networks))
@info extrema(links.(networks))

#=
!!! warning "Paging and API queries"

    Queries using the *mangal* API use a paging system. This dataset only has 50
    networks, which happens to fit in a single page, but this may not be the case
    for all datasets. Go check out the appropriate documentation.
=#

# Note that this step is constrained in part by queries over the network, as we
# are requesting about 50 different networks. In order to simplify the analysis,
# we will make these networks binary, using [`render`](@ref):

networks = [render(Binary, network) for network in networks]

@info extrema(richness.(networks))
@info extrema(links.(networks))

# To simplify the analysis, we will focus on a single property of interest (the
# spectral radius), and look at the distribution of the value for 99 replicates
# -- this is not enough for an actualy analysis, but good enough for a
# demonstration.

function nm_comparison(N::SpeciesInteractionNetwork)
    R = [structuralmodel(NicheModel, N) for _ in 1:99]
    return spectralradius.(R)
end

# We will simplify things further by onlly looking at the *z*-score, so we can
# extract the mean and standard deviation. We will return a *function* to
# calculate the z-score instead of the raw values:

summarizer(x) = (x₀) -> (x₀ - Statistics.mean(x)) / Statistics.std(x)

# And the *z*-scores are:

Z = zeros(Float64, length(networks))
for (i,network) in enumerate(networks)
    zₙ = summarizer(nm_comparison(network))
    Z[i] = zₙ(spectralradius(network))
end

# How well does the niche model approximates the spectral radius? We can look at
# the distribution of z-scores, and their relationship with *e.g.* connectance:

f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax1 = CairoMakie.Axis(f[1,1], xlabel="Z-score", ylabel="Probability")
ax2 = CairoMakie.Axis(f[1,2], xlabel="Connectance", ylabel="Z-score")
CairoMakie.hist!(ax1, Z; normalization=:probability, fillto=0.0, color=(:slategray, 0.4), bins=20)
CairoMakie.scatter!(ax2, connectance.(networks), Z; color=:slategray)
CairoMakie.vlines!(ax1, [0.0]; color=:black)
CairoMakie.hlines!(ax2, [0.0]; color=:black)
CairoMakie.tightlimits!(ax1)
CairoMakie.tightlimits!(ax2)
CairoMakie.current_figure()

# It seems that the niche model tends to under-estimate the spectral radius,
# although this is more marked for less densely connected networks.