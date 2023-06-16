# Beta-diversity

!!! abstract

    These methods are used to measure the beta-diversity of two networks, by partitioning variation into the dissimilarity of species (βS), interactions (βOS), and whole networks (βWN).

## Components of networks beta diversity

```@docs
BetaDivComponent
βS
βOS
βWN
```

## Beta diversity measures

The [`betadiversity`](@ref) function will *always* return a named tuple with three entries, named

```@docs
betadiversity
```