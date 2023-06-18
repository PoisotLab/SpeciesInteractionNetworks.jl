# Nestedness

!!! abstract

    The methods presented in this page measure the nestedness of a network. Nestedness is usually restricted to biparite networks, although following the arguments laid out by [Staniczenko2013ghost](@citet), we consider [`spectralradius`](@ref) to be a measure of nestedness.

## η

!!! warning "Degree distribution and η"

    The η measure of nestedness is *invariant* for a given degree distribution. In other words, two networks with the same (joint) degree distribution will *always* have the same value of η. As a result, network permutations using [`swap!`](@ref) and a constraint on the degree will not generate an appropriate null sample. This is also true when only one side on the degree distribution is maintained, for the measure of nestedness on this side.

```@docs
η
```

## NODF

## Spectral radius

```@docs
spectralradius
```
