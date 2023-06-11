# Permutations of networks

!!! abstract

    The methods presented in this page perform network *permutation*, *i.e.* they move interactions around while also respecting a number of constraints. Permutations are used in null hypothesis testing, or can be used alongside *e.g.* simulated annealing to bring networks closer to a specified structure.

The functions for permutations are using an edge-swap algorithm, in which the
endpoint of interactions is switched to re-wire the network without changing the degree distribution. Each call to the `swap!` function will *modify* the network, and perform a single edge swap.

!!! info "Maximum iterations when doing permutations"

    Every permutation will try up to `SpeciesInteractionNetworks.SWAP_MAXITER` times (defaults to 100) to find a suitable pair of edges to swap, and then return the network *unshuffled* if they failed to find a suitable pair of edges to swap. This value *can* be changed.

## Permutation constraints

Permutations are constrained, in that we can guarantee that the resulting network may have structural properties that are similar to the original network. The type of constraint we apply is determined by the `PermutationConstrant` enumerated type.

```@docs
PermutationConstraint
```

## Permutation of a network

Note that the permutations are currently limited to networks with `Binary` interactions.
