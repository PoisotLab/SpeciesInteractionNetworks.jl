# # Measuring spcificity

# !!! example "Measuring specificity"
#     
#     In this example, we will...

# We can generate an example network with three different degrees of specificity:

using SpeciesInteractionNetworks

nodes = Bipartite([:A, :B, :C, :D, :E], [:a, :b, :c, :d, :e, :f])
edges = Quantitative([1 0 0 0 0; 2 0 0 0 0; 1 1 1 0 0; 4 3 2 1 0; 4 4 4 3 0])
N = SpeciesInteractionNetwork(nodes, edges)

# We can calculate the specificity of the top-level species:

spe_scores = specificity(N)

# The output of specificity is a dictionary, where the `species(N, 1)` are keys,
# and the score for each of these species are the values. We can, for example,
# look at the specificity for species `:D`:

spe_scores[:D]

# Note that if we want the value for a smaller number of species, it is faster
# to call the function with a single species name:

specificity(N, :D)

# !!! info "Making sense of the score"
# 
#     The Paired Differences Index will *always* return values in the unit interval, and these values are independent from one species to the next. In the example above, species `:A` and `:B` have the same (maximum) specificity because they use a single resource. The purpose of the Paired Differences Index is to express specificity in a way that is *not* affected by the total interaction strenght of the species, because what is well understood can be measured without confounders.

# As always, keep in mind that the ordering of keys in the dictionary is not
# fixed. Therefore, it is probably safer to iterate over the `species(N, 1)`
# when looking for specific values.
