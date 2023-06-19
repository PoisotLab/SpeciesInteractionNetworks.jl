# Random draws

!!! abstract

    This page presents an overview of methods to draw from a probabilistic network, and ways to generate structural null models from a binary network under various constraints.

## Drawing from a probabilistic network

```@docs
randomdraws
```

## Null models for hypothesis testing

The usual null models of network structure can be generated using a
[`PermutationConstraint`](@ref) as the first argument. Although technically
speaking, these are not permutations, it is nevertheless useful to map the
constraints from one method of network generation to another.

```@docs
nullmodel
```

## Species-level contributions based on null models

The null model approach can be extended to the contribution of each species. We
offer a generic interface to generating pseudo-random networks where the target
species has its interactions randomized, while all other interactions remain
fixed.

```@docs
speciescontribution
```

## Linear filtering

The linear filter from [Stock2017Linear](@citet) is a general formulation of the
null models presented above.

```@docs
linearfilter
```
