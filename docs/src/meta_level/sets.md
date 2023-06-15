# Networks as sets

!!! abstract

    TODO

## Operations

```@docs
Base.union(U::T, V::T) where T<:(SpeciesInteractionNetwork{<:SpeciesInteractionNetworks.Partiteness, <:Binary})
```

```@docs
Base.intersect(::SpeciesInteractionNetwork,::SpeciesInteractionNetwork)
Base.setdiff(::SpeciesInteractionNetwork,::SpeciesInteractionNetwork)
```