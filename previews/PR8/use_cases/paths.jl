# # Shortest paths

# !!! example "Shortest paths"
#     
#     In this example, we will look at the simplified food web of the purple pitcher plant *Sarracenia purpurea*, as an excuse to see how the functions to handle shortest path analyses work.

# The food web of the purple pitcher plan is surprisingly simple. In this
# demonstration, we will use the version presented in
# [Ellison2021Regulation](@citet), itself an adaptation from
# [Baiser2013Predicting](@citet).

using SpeciesInteractionNetworks

# We will create the nodes first:

nodes = Unipartite([
    :Fletcherimyia,
    :Wyeomyia,
    :Protozoa,
    :Habrotrocha,
    :Bacteria,
    :Sarraceniopus,
    :POM,
    :Metriocnemus,
    :Detritus
])

# In order to facilitate our work, we will start with an empty food web, and
# fill in interactions later:

edges = Binary(zeros(Bool, (richness(nodes), richness(nodes))))
N = SpeciesInteractionNetwork(nodes, edges)

# The interactions can be added one by one. Note that we represent interactions
# using the *from* → *to* syntax.

N[:Fletcherimyia, :Fletcherimyia] = true
N[:Fletcherimyia, :Protozoa] = true
N[:Fletcherimyia, :Habrotrocha] = true
N[:Fletcherimyia, :Wyeomyia] = true
N[:Fletcherimyia, :Bacteria] = true
N[:Fletcherimyia, :Detritus] = true
N[:Wyeomyia, :Protozoa] = true
N[:Wyeomyia, :Habrotrocha] = true
N[:Wyeomyia, :Bacteria] = true
N[:Habrotrocha, :Protozoa] = true
N[:Habrotrocha, :Bacteria] = true
N[:Habrotrocha, :POM] = true
N[:Protozoa, :Bacteria] = true
N[:Bacteria, :POM] = true
N[:Bacteria, :Detritus] = true
N[:Sarraceniopus, :Bacteria] = true
N[:Sarraceniopus, :POM] = true
N[:POM, :Metriocnemus] = true # This one is strange but let's roll with it
N[:Metriocnemus, :Detritus] = true

N

# With this network, we can start by measuring the length of the shortest path
# between Fletcherimyia and all other nodes of the food web:

shortestpath(N, :Fletcherimyia)

# By default, the [`shortespath`](@ref) function will use the
# [`BellmanFord`](@ref) algorithm. It is theoretically slower than
# [`Dijkstra`](@ref), but in practice, we use an early termination check which
# makes it compare very favorably.

# We can get the path between Fletcherimyia and, for example, Metriocnemus,
# using [`pathbetween`](@ref):

pathbetween(N, :Fletcherimyia, :Metriocnemus)

# This returns a list of interactions. Note that the [`pathbetween`](@ref)
# method will return a *single* path, even though several may exist.
