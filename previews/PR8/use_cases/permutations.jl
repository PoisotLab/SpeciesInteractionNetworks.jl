# # Permutations

# To showcase [`swap!`](@ref) in practice, we will work through through a simple
# example of (i) generating a perfectly nested network, (ii) shuffling
# interactions by maintaining the generality of top-level species, and (iii)
# looking at the way the nestdeness of the entire network changes with each
# successive swap.

using SpeciesInteractionNetworks
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide

# We can generate a nested network rather easily, by creating a matrix of binary
# interactions, where the species interact with species from a lower rank:

A = zeros(Bool, (10, 14))
for i in axes(A, 1)
    for j in axes(A, 2)
        if i <= j
            A[i,j] = true
        end
    end
end

# We can declare a network *without* having to define all of the species, by
# first wrapping our matrix inside a `Binary` type, and then generating a
# `Bipartite` species set with the right number of species:

edges = Binary(A)
nodes = Bipartite(edges)
N = SpeciesInteractionNetwork(nodes, edges)

# The initial nestedness of this network, measured using [η](@ref) is (network,
# top-level contribution, bottom-level contribution):

(η(N), η(N,1), η(N,2))

# In order to generate the series of successive permutations, we will define an
# empty array of values, and then for each successive step, calculate the
# nestedness of the network, and then swap interactions under the given
# constraint.

nestedness_series = zeros(Float64, 1000)
for i in axes(nestedness_series, 1)
    nestedness_series[i] = η(N)
    swap!(N, Generality)
end

#!!! warning "A note about swaps and underlying nodes/edges"
#
#    When we perform the `swap!` operation, we are modifying the network (this is
#    what we want!), but we are also modifying the `edges` object. If you want to
#    re-use the edges in another network, be mindful of the fact that this will
#    be the *randomized* edges. See [`copy`](@ref) for a way to create new copies of a network.

# Finally, we can plot the result, to check that 1000 swaps are enough to bring
# us to some sort of equilibrium of the randomized nestedness:

f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Swap", ylabel = "Nestedness")
CairoMakie.lines!(ax, nestedness_series, color=(:black, 0.5))
CairoMakie.tightlimits!(ax)
CairoMakie.current_figure()
