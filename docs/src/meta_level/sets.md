# Networks as sets

!!! abstract

    A number of operations can be applied to networks. This include `union`, `intersect`, and `setdiff`. These are useful when workin on network beta-diversity.

## Operations

```@docs
Base.union(U::T, V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
Base.intersect(U::T, V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
Base.setdiff(U::T, V::T) where {T <: SpeciesInteractionNetwork{<:Partiteness, <:Binary}}
```
