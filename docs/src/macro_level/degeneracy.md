# Disconnected species

!!! abstract

    When generating networks, in particular when using [`randomdraws`](@ref) following [`nullmodel`](@ref), there is a chance that some species will have no interactions. In some cases, it may be relevant to identify and remove these species. The methods presented in this page offer a way to do this.

## Detection of disconnected species

```@docs
isdisconnected
isdegenerate
```

## Simplification of networks with disconnected species

```@docs
simplify
```

!!! question "Where is `simplify!`"

    In the `EcologicalNetworks.jl` package, the `simplify!` function would drop species from a network. In this package, because we try to keep networks as immutable as possible, we have not implemented this function. If you really want to reproduce this behavior, you can wrap whatever command will create the network in `simplify`.
