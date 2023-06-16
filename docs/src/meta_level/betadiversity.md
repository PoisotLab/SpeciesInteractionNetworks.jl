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

!!! tip "What about the species turnover component?"

    In [Poisot2012dissimilarity](@cite), we introduced the idea that the impact of species turnover can often be expressed as the difference between the whole-network ([`βWN`](@ref)) and overlapping-species ([`βOS`](@ref)) dissimilarities. This is only true for *some* measures. Following the arguments laid out in [Poisot2022Dissimilarity](@cite), we have *not* added this as a built-in function. If there is a need to measure the impact of turnover, it is recommended to express it as (wn-os)/wn. 

## Beta diversity measures

The [`betadiversity`](@ref) function will *always* return a named tuple with three entries, named

```@docs
betadiversity
```