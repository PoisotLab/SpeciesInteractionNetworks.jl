# Measures of centrality

!!! abstract

    Centrality can help in quantifying the importance of species in a network. These function will measure the centrality of all species under different algorithms. There is a single wrapper function called `centrality`, which uses an optional first argument to specify the algorithm to use.

The centrality scores are returned so that they *sum* to one. This is intended
to make sure that *within a network*, the values for different nodes are
comparable.

## Implemented algorithms

```@docs
CentralityMethod
```

```@docs
KatzCentrality
EigenvectorCentrality
ClosenessCentrality
ResidualClosenessCentrality
GeneralizedClosenessCentrality
```

## Centrality function

```@docs
centrality
```

## Illustration

```@example 1
using SpeciesInteractionNetworks

nodes = Unipartite([:Wyeomyia, :Metriocnemus, :Fletcherimyia, :Habrotrocha, :Protozoa, :Sarraceniopus, :Bacteria, :Insects])
edges = Binary(zeros(Bool, (richness(nodes), richness(nodes))))
N = SpeciesInteractionNetwork(nodes, edges)

N[:Metriocnemus, :Insects] = true
N[:Fletcherimyia, :Insects] = true
N[:Bacteria, :Insects] = true
N[:Sarraceniopus, :Bacteria] = true
N[:Habrotrocha, :Bacteria] = true
N[:Protozoa, :Bacteria] = true
N[:Wyeomyia, :Habrotrocha] = true
N[:Wyeomyia, :Protozoa] = true
N[:Wyeomyia, :Bacteria] = true

interactions(N)
```

change alpha

```@example 1
attenuation = 10.0.^LinRange(-3, 0, 10)
c_insect = zeros(length(attenuation))
c_bacteria = zeros(length(attenuation))
c_protozoa = zeros(length(attenuation))
for (i,α) in enumerate(attenuation)
    ci = centrality(KatzCentrality, N; α=α)
    c_insect[i] = ci[:Insects]
    c_bacteria[i] = ci[:Bacteria]
    c_protozoa[i] = ci[:Protozoa]
end
```

plot

```@example 1
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide

f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Attenuation", ylabel = "Centrality", xscale=log10)
CairoMakie.lines!(ax, attenuations, c_insect, color=(:black, 0.5), label="Insect")
CairoMakie.lines!(ax, attenuations, c_bacteria, color=(:green, 0.5), label="Bacteria")
CairoMakie.lines!(ax, attenuations, c_protozoa, color=(:orange, 0.5), label="Protozoa")
CairoMakie.tightlimits!(ax)
CairoMakie.axislegend()
CairoMakie.current_figure()
```