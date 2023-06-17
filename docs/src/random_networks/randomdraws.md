# Random draws

!!! abstract

    This page presents an overview of methods to draw from a probabilistic network, and ways to generate structural null models from a binary network under various constraints. We present an illustration using plant-pollinator data [Dupont2003Structure](@cite) to show how these functions can be applied to test the significance of network structure.

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

## Linear filtering

The linear filter from [Stock2017Linear](@citet) is a general formulation of the
null models presented above.

```@docs
linearfilter
```

## Illustration

In order to illustrate the use of null models, we will look at the data from
[Dupont2003Structure](@citet), and specifically generate multiple networks under
different null models, then compare the nestedness (using [`η`](@ref)) to the
empirical network.

```@example 1
using SpeciesInteractionNetworks
import DelimitedFiles
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide
import Statistics
```

The data are available from the
[IWDB](http://www.ecologia.ib.usp.br/iwdb/html/dupont_et_al.html), but in order
to avoid making unecessary calls to their webserver, we have reproduced a
version here:

```@example 1
int_mat = Bool[
    1 1 1 0 1 1 0 0 1 0 0; 1 0 0 1 1 1 0 1 0 1 0; 1 0 1 0 1 1 1 0 0 0 0;
    0 1 0 1 0 1 1 1 0 1 0; 1 0 1 0 1 1 1 0 0 1 0; 0 1 1 0 1 0 1 0 1 0 1;
    1 1 1 1 0 0 1 0 0 0 0; 1 0 0 0 1 1 0 0 0 1 0; 1 1 1 0 0 0 0 1 0 0 0;
    1 1 0 0 0 1 1 0 0 0 0; 1 1 0 1 0 0 0 0 1 0 0; 1 1 1 0 0 0 0 1 0 0 0;
    0 1 0 1 0 0 0 1 0 0 0; 1 0 1 0 1 0 0 0 0 0 0; 0 0 1 0 1 0 0 0 0 0 1;
    1 0 1 1 0 0 0 0 0 0 0; 0 0 1 1 1 0 0 0 0 0 0; 0 1 0 0 0 0 1 0 1 0 0;
    0 0 0 0 0 1 1 0 0 0 0; 0 0 0 1 0 1 0 0 0 0 0; 1 0 0 0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 1 0 0 0; 1 0 0 0 0 0 0 0 1 0 0; 0 0 1 0 0 0 0 0 1 0 0;
    0 0 0 1 0 0 0 0 0 0 0; 0 1 1 0 0 0 0 0 0 0 0; 1 0 0 1 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0
]
nothing # hide
```

We can turn this into a network (without species names!):

```@example 1
edges = Binary(int_mat)
nodes = Bipartite(edges)

N = SpeciesInteractionNetwork(nodes, edges)
@info "$(richness(N,1)) pollinators"
@info "$(richness(N,2)) plants"
@info "$(length(N)) interactions"
```


We first measure the nestedness of the network:

```@example 1
n0 = η(N)
```

The next step is to generate a template probabilistic network under a specific
null model. Here, we will focus on the null model based on the joint degree
distribution, as used by *e.g.* [Bascompte2003nested](@citet):

```@example 1
Nd = nullmodel(Degree, N)
```

We can draw samples from this network, and measure their nestedness:

```@example 1
Rd = [randomdraws(Nd) for _ in 1:999]
nd = η.(Rd)
@info round.((minimum(nd), Statistics.median(nd), maximum(nd)); digits=4)
```

This is all we need to plot the results:

```@example 1
f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Nestedness", ylabel="Probability")
CairoMakie.hist!(ax, nd; normalization=:probability, fillto=0.0, color=(:slategray, 0.4), bins=20)
CairoMakie.vlines!(ax, [n0], color=:black, linestyle=:dash)
CairoMakie.tightlimits!(ax)
CairoMakie.ylims!(ax, (0.0, 0.2))
CairoMakie.current_figure()
```

In practice, we are often interested in deriving a *p*-value from the comparison
of the empirical and null values of the structure measure. Note that,
functionally, the generation of null models can be seen as permutation testing,
and therefore we can approximate the *p*-value corresponding to the hypothesis
that the network is *more* nested than expected under its degree distribution by
looking at the proportion of randomized values that are larger than the
empirical observation:

```@example 1
@info "𝑝 ≈ $(round(Statistics.mean(n0 .<= nd); digits=3))"
```
