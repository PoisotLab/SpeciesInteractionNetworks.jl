# Permutations of networks

!!! abstract

    The methods presented in this page perform network *permutation*, *i.e.* they move interactions around while also respecting a number of constraints. Permutations are used in null hypothesis testing, or can be used alongside *e.g.* simulated annealing to bring networks closer to a specified structure.

The functions for permutations are using an edge-swap algorithm, in which the
endpoint of interactions is switched to re-wire the network without changing the
degree distribution. Each call to the `swap!` function will *modify* the
network, and perform a single edge swap.

!!! info "Maximum iterations when doing permutations"

    Every permutation will try up to `SpeciesInteractionNetworks.SWAP_MAXITER` times (defaults to 100) to find a suitable pair of edges to swap, and then return the network *unshuffled* if they failed to find a suitable pair of edges to swap. This value *can* be changed.

## Permutation constraints

Permutations are constrained, in that we can guarantee that the resulting
network may have structural properties that are similar to the original network.
The type of constraint we apply is determined by the `PermutationConstrant`
enumerated type.

```@docs
PermutationConstraint
```

## Permutation of a network

Note that the permutations are currently limited to networks with `Binary`
interactions.

```@docs
swap!
```

##  Illustration

To showcase `swap!` in practice, we will work through through a simple example
of (i) generating a perfectly nested network, (ii) shuffling interactions by
maintaining the generality of top-level species, and (iii) looking at the way
the nestdeness of the entire network changes with each successive swap.

```@example 1
using SpeciesInteractionNetworks
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide
```

We can generate a nested network rather easily, by creating a matrix of binary
interactions, where the species interact with species from a lower rank:

```@example 1
A = zeros(Bool, (10, 14))
for i in axes(A, 1)
    for j in axes(A, 2)
        if i <= j
            A[i,j] = true
        end
    end
end
```

We can declare a network *without* having to define all of the species, by first
wrapping our matrix inside a `Binary` type, and then generating a `Bipartite`
species set with the right number of species:

```@example 1
edges = Binary(A)
nodes = Bipartite(edges)
N = SpeciesInteractionNetwork(nodes, edges)
```

The initial nestedness of this network, measured using [η](@ref) is (network,
top-level contribution, bottom-level contribution):

```@example 1
(η(N), η(N,1), η(N,2))
```

In order to generate the series of successive permutations, we will define an
empty array of values, and then for each successive step, calculate the
nestedness of the network, and then swap interactions under the given
constraint.

```@example 1
nestedness_series = zeros(Float64, 1000)
for i in axes(nestedness_series, 1)
    nestedness_series[i] = η(N)
    swap!(N, Generality)
end
```

!!! warning "A note about swaps and underlying nodes/edges"

    When we perform the `swap!` operation, we are modifying the network (this is
    what we want!), but we are also modifying the `edges` object. If you want to
    re-use the edges in another network, be mindful of the fact that this will
    be the *randomized* edges. See [`copy`](@ref) for a way to create new copies of a network.

Finally, we can plot the result, to check that 1000 swaps are enough to bring us
to some sort of equilibrium of the randomized nestedness:

```@example 1
f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Swap", ylabel = "Nestedness")
CairoMakie.lines!(ax, nestedness_series, color=(:black, 0.5))
CairoMakie.tightlimits!(ax)
CairoMakie.current_figure()
```
