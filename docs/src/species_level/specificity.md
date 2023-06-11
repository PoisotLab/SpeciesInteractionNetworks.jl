# Specificity

!!! abstract

    The specificity of species in a network is measured either as a function of the proportion of resouces they effectively use (for binary networks), or as a function of the distribution of their performance on these resources (for quantitative networks).

## Measures of specificity

The packages relies on the *Paired Difference Index* to calculate specificity,
as it can be applied to both binary and quantitative data.

```@docs
specificity
```

## Illustration

We can generate an example network with three different degrees of specificity:

```@example 1
using SpeciesInteractionNetwork

nodes = Bipartite([:A, :B, :C, :D], [:a, :b, :c, :d, :e, :f])
edges = Quantitative([1 0 0 0 0; 2 0 0 0 0; 1 1 1 0 0; 4 3 2 1 0])
N = SpeciesInteractionNetwork(nodes, edges)
```

We can calculate the specificity of the top-level species:

```@example 1
spe_scores = specificity(N)
```

!!! info "Making sense of the score"

    The Paired Differences Index will *always* return values in the unit interval, and these values are independent from one species to the next. In the example above, species `:A` and `:B` have the same (maximum) specificity because they use a single resource. The purpose of the Paired Differences Index is to express specificity in a way that is *not* affected by the total interaction strenght of the species, because what is well understood can be measured without confounders.

The output of specificity is a dictionary, where the `species(N, 1)` are keys,
and the score for each of these species are the values. We can, for example,
look at the specificity for species `:D`:

```@example 1
spe_scores[:D]
```