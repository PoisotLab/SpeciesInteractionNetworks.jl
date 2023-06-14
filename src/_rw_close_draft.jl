using SpeciesInteractionNetworks

edges = Probabilistic(rand(Float64, (20, 20)))
nodes = Unipartite(edges)
N = SpeciesInteractionNetwork(nodes, edges)

A = Matrix(N.edges.edges)
m = A ./ sum(A; dims=2)

using LinearAlgebra

Hdotj = (I - m[1:end.!=j,1:end.!=j])^(-1) * ones(size(m,1)-1)

