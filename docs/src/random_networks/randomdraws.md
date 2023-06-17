# Random draws

!!! abstract

    The methods in this page draw random binary networks from probabilistic networks

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
import Downloads
import DelimitedFiles
import CairoMakie
import Statistics

int_mat_path = Downloads.download("http://www.ecologia.ib.usp.br/iwdb/data/plant_pollinator/text/dupont_matr.txt")
int_mat = Bool.(DelimitedFiles.readdlm(int_mat_path))

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
f = CairoMakie.Figure()
ax = CairoMakie.Axis(f[1,1], xlabel="Nestedness", ylabel="Samples")
CairoMakie.hist!(ax, nd; normalization=:probability, fillto=0.0, color=(:blue, 0.4))
CairoMakie.vlines!(ax, [n0], color=:black)
CairoMakie.tightlimits!(ax)
CairoMakie.current_figure()
```
