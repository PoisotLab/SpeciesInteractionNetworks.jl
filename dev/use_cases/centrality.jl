# # Measuring centrality

# !!! example "Centrality"
#     
#     In this example, we will...

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
N[:Protozoa, :Habrotrocha] = true
N[:Wyeomyia, :Habrotrocha] = true
N[:Wyeomyia, :Protozoa] = true
N[:Wyeomyia, :Bacteria] = true

interactions(N)

# Change α

attenuation = 10.0.^LinRange(-1, 0, 20)
c_insect = zeros(length(attenuation))
c_bacteria = zeros(length(attenuation))
c_protozoa = zeros(length(attenuation))
for (i,α) in enumerate(attenuation)
    ci = centrality(KatzCentrality, N; α=α)
    c_insect[i] = ci[:Insects]
    c_bacteria[i] = ci[:Bacteria]
    c_protozoa[i] = ci[:Protozoa]
end

# Make a plot

import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide

f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Attenuation", ylabel = "Centrality", xscale=log10)
CairoMakie.lines!(ax, attenuation, c_insect, color=(:black, 0.5), label="Insect")
CairoMakie.lines!(ax, attenuation, c_bacteria, color=(:green, 0.5), label="Bacteria")
CairoMakie.lines!(ax, attenuation, c_protozoa, color=(:orange, 0.5), label="Protozoa")
CairoMakie.tightlimits!(ax)
CairoMakie.axislegend(position=:lb)
CairoMakie.current_figure()


# same with closeness centrality

attenuation = 10.0.^LinRange(-1, 0, 20)
c_insect = zeros(length(attenuation))
c_bacteria = zeros(length(attenuation))
c_protozoa = zeros(length(attenuation))
for (i,α) in enumerate(attenuation)
    ci = centrality(GeneralizedClosenessCentrality, N; α=α)
    c_insect[i] = ci[:Insects]
    c_bacteria[i] = ci[:Bacteria]
    c_protozoa[i] = ci[:Protozoa]
end

# Make a plot

import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide

f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Attenuation", ylabel = "Centrality", xscale=log10)
CairoMakie.lines!(ax, attenuation, c_insect, color=(:black, 0.5), label="Insect")
CairoMakie.lines!(ax, attenuation, c_bacteria, color=(:green, 0.5), label="Bacteria")
CairoMakie.lines!(ax, attenuation, c_protozoa, color=(:orange, 0.5), label="Protozoa")
CairoMakie.tightlimits!(ax)
CairoMakie.axislegend(position=:lb)
CairoMakie.current_figure()