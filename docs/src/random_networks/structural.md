# Structural models

!!! abstract

    These functions generate random networks based on structural models. The models are given as types (specifically subtypes of [`StructuralModel`](@ref)), and return a network the type of which depends on the generative algorithm, when passed to the [`structuralmodel`](@ref) function.

## Models

```@docs
StructuralModel
NicheModel
```

## The generative function

```@docs
structuralmodel
```

## Illustration

In this example, we will look at the relationship between connectance and
spectral radius ([`spectralradius`](@ref)) for randomly generated networks under
the niche model.

```@example 1
using SpeciesInteractionNetworks
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide
```

We can generate a number of random networks:

```@example 1
R = [structuralmodel(NicheModel, 25, rand()*0.5) for _ in 1:999]
nothing # hide
```

We can get the connectance of the random networks -- recall that the niche model
is not going to strictly respect the target connectance:

```@example 1
co = connectance.(R)
nothing # hide
```

Finally, we can look at the expected relationship between connectance and
spectral radius:

```@example 1
f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Connectance", ylabel="Spectral radius")
CairoMakie.scatter!(ax, co, spectralradius.(R; correction=:connectance), color=(:slategray,0.4))
CairoMakie.tightlimits!(ax)
CairoMakie.xlims!(ax, (0.0, 0.5))
CairoMakie.current_figure()
```
