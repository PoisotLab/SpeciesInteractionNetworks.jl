# # Spectral radius of the niche model

# !!! example "The niche model of food webs"
#     
#     In this example, we will...

using SpeciesInteractionNetworks
import CairoMakie
CairoMakie.activate!(px_per_unit=2) #hide

# We can generate a number of random networks:

R = [structuralmodel(NicheModel, 25, rand()*0.5) for _ in 1:999]

# We can get the connectance of the random networks -- recall that the niche model
# is not going to strictly respect the target connectance:

co = connectance.(R)

# Finally, we can look at the expected relationship between connectance and
# spectral radius:

f = CairoMakie.Figure(backgroundcolor = :transparent, resolution = (800, 300))
ax = CairoMakie.Axis(f[1,1], xlabel="Connectance", ylabel="Spectral radius")
CairoMakie.scatter!(ax, co, spectralradius.(R), color=(:slategray,0.4))
CairoMakie.tightlimits!(ax)
CairoMakie.xlims!(ax, (0.0, 0.5))
CairoMakie.current_figure()
