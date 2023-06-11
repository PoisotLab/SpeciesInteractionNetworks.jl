"""
    CentralityMethod

All algorithms for centrality are subtypes of `CentralityMethod`.
"""
abstract type CentralityMethod end

abstract type ClosenessCentrality <: CentralityMethod end
abstract type DegreeCentrality <: CentralityMethod end
abstract type HarmonicCentrality <: CentralityMethod end
abstract type EigenvectorCentrality <: CentralityMethod end

"""
    KatzCentrality

This type is used to perform the Katz centrality analysis.
"""
abstract type KatzCentrality <: CentralityMethod end

"""
    centrality(::Type{KatzCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary, Probabilistic}}; α::AbstractFloat=0.1)

This measure gives a different weight to every subsequent connection (`α`). `α`
is a weight, specifically the weight of each subsequent move away from the node,
and therefore must be positive.

Note that internally, this function uses linear algebra shortcuts (rather than a
sum over path lengths) to calculate the centrality, which is faster.

The values are returned as a dictionary, and the values are futher ranged so
that the sum of centralities is equal to one.

###### References

[Katz1953new](@cite)

[Junker2008Analysis](@cite)
"""
function centrality(::Type{KatzCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary, Probabilistic}}; α::AbstractFloat=0.1)    
    @assert 0.0 <= α <= 1.0

    A = Matrix(N.edges.edges)
    C = (inv(I - α*A')-I)*ones(richness(N))

	return Dict(zip(species(N), C ./ sum(C)))
end

@testitem "We cannot use Katz centrality with a outside the unit interval" begin
    nodes = Unipartite([:A, :B, :C])
    edges = Binary(Bool[0 1 1; 1 0 1; 0 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test_throws AssertionError centrality(KatzCentrality, N; α = 1.2)
    @test_throws AssertionError centrality(KatzCentrality, N; α = -0.5)
end

#=

"""
    centrality_closeness(N::UnipartiteNetwork; nmax::Int64=20)

The function calls `shortest_path` internally -- the `nmax` argument is the
maximal path length that will be tried.

#### References

- Bavelas, A., 1950. Communication Patterns in Task‐Oriented Groups. The Journal
  of the Acoustical Society of America 22, 725–730. doi:10.1121/1.1906679
"""
function centrality_closeness(N::UnipartiteNetwork; nmax::Int64=20)
  d = shortest_path(N, nmax=nmax)
  n = richness(N)-1
  d[diagind(d)] .= 0
  interm = sum(d; dims=2)
  interm = vec(n ./ interm)
  for i in eachindex(interm)
    interm[i] = interm[i] == Inf ? 0.0 : interm[i]
  end
  return Dict(zip(species(N), interm))
end

"""
    centrality_harmonic(N::UnipartiteNetwork; nmax::Int64=20)

The function calls `shortest_path` internally -- the `nmax` argument is the
maximal path length that will be tried.
"""
function centrality_harmonic(N::UnipartiteNetwork; nmax::Int64=20)
  d = shortest_path(N, nmax=nmax)
  n = richness(N)-1
  d[diagind(d)] .= 0
  interm = 1.0 ./ d
  for i in eachindex(interm)
    interm[i] = isinf(interm[i]) ? 0.0 : interm[i]
  end
  interm = sum(d; dims=2)./n
  return Dict(zip(species(N), interm))
end

"""
    centrality_degree(N::UnipartiteNetwork)

Degree centrality, corrected by the maximum degree (the most central species has
a degree of 1).
"""
function centrality_degree(N::UnipartiteNetwork)
  d = degree(N)
  dm = maximum(values(d))
  return Dict([p.first=>p.second/dm for p in d])
end


"""
    centrality_eigenvector(N::AbstractUnipartiteNetwork)

Eigen centrality, corrected by the maximum degree (the most central species has
a degree of 1).
"""
function centrality_eigenvector(N::AbstractUnipartiteNetwork)
  evals, evecs = eigen(Array(N.edges))
  emax, epos = findmax(real.(evals))
  cvals = vec(real.(evecs[:,epos]))
  cvals ./= maximum(cvals)
  return Dict(zip(species(N), cvals))
end
=#