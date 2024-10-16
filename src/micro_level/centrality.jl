"""
    CentralityMethod

All algorithms for centrality are subtypes of `CentralityMethod`. These
algorithms do not take additional arguments, which are instead passed to the
[`centrality`](@ref) method.
"""
abstract type CentralityMethod end

"""
    ClosenessCentrality

Closeness centrality is defined as the reciprocal of farness, *i.e.* 1 over the
sum of shortest paths from the source node to all other nodes. In the case of
probabilistic networks, we use the random walk approximation of closeness
centrality, where the distance is defined as first passage time.

###### References

[Bavelas1950Communication](@citet*)
"""
abstract type ClosenessCentrality <: CentralityMethod end

"""
    ResidualClosenessCentrality

###### References

[Dangalchev2006Residual](@citet*)
"""
abstract type ResidualClosenessCentrality <: CentralityMethod end

"""
    GeneralizedClosenessCentrality

###### References

[Dangalchev2011Residual](@citet*)
"""
abstract type GeneralizedClosenessCentrality <: CentralityMethod end

"""
    KatzCentrality

This type is used to perform the Katz centrality analysis.

###### References

[Katz1953new](@citet*)
"""
abstract type KatzCentrality <: CentralityMethod end

"""
    EigenvectorCentrality

This type is used to perform the Eigenvector centrality analysis.
"""
abstract type EigenvectorCentrality <: CentralityMethod end

"""
    centrality(N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary, Probabilistic}})

If no type is given as the first argument, the centrality function will use
[`ClosenessCentrality`](@ref).
"""
function centrality(N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary, Probabilistic}})
    return centrality(ClosenessCentrality, N)
end

"""
    centrality(::Type{ClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Probabilistic})

The closeness centrality of a probabilistic network is measured using the
*random walk* closeness centrality, where the distance between two nodes i and j
is measured as the first passage time on node j of a random walk starting on
node i. This function does not look for shortest paths, and is therefore
reasonably fast.
"""
function centrality(::Type{ClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Probabilistic})
    A = Matrix(N.edges.edges)
    m = A ./ sum(A; dims=2)
    e = ones(size(m,1)-1)
    ci = zeros(richness(N))
    for i in 1:richness(N)
        Hdoti = (I - m[1:end.!=i,1:end.!=i])^(-1) * e
        ci[i] = n/sum(Hdoti)
    end
    return Dict(zip(species(N), ci ./ sum(ci)))
end

"""
    centrality(::Type{ClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

The closeness centrality of a binary network is the reciprocal of the sum of
shortest paths.
"""
function centrality(::Type{ClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    ci = zeros(richness(N))
    for i in 1:richness(N)
        sp = shortestpath(BellmanFord, N, species(N)[i])
        if isempty(sp)
            ci[i] = 0.0
        else
            ci[i] = 1.0 / sum(collect(values(sp))) * (richness(N)-1)
        end
    end
    return Dict(zip(species(N), ci ./ sum(ci)))
end

"""
    centrality(::Type{ResidualClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})

The residual closeness centrality of a binary network is the sum of the
residuals of two to the power of the length of shortest paths.
"""
function centrality(::Type{ResidualClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary})
    ci = zeros(richness(N))
    for i in 1:richness(N)
        sp = shortestpath(BellmanFord, N, species(N)[i])
        if isempty(sp)
            ci[i] = 0.0
        else
            ci[i] = sum(1.0/(2.0.^values(sp)))
        end
    end
    return Dict(zip(species(N), ci ./ sum(ci)))
end

"""
    centrality(::Type{GeneralizedClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary}; α=0.1)

The generalized closeness centrality of a binary network ranges from the degree
centrality (α is 0) to the number of reachable nodes (α=1).
"""
function centrality(::Type{GeneralizedClosenessCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Binary}; α=0.1)
    ci = zeros(richness(N))
    for i in 1:richness(N)
        sp = shortestpath(BellmanFord, N, species(N)[i])
        if isempty(sp)
            ci[i] = 0.0
        else
            ci[i] = sum(α.^values(sp))
        end
    end
    return Dict(zip(species(N), ci ./ sum(ci)))
end

"""
    centrality(::Type{KatzCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary, Probabilistic}}; α::AbstractFloat=0.1)

This measure gives a different weight to every subsequent connection (`α`). `α`
is a weight, specifically the attenuation of each subsequent move away from the
node, and therefore must be positive.

Note that internally, this function uses linear algebra shortcuts (rather than a
sum over path lengths) to calculate the centrality, which is faster. As a
consequence, there is, for each network, a maximal value of α that is the
reciprocal of the absolute value of the leading eigenvalue of the adjacency
matrix. This α has no analytical solution when the adjacency matrix has complex
eigenvalues. It is recommended to use this centrality algorithm as part of a
`try`/`catch` block if using large values of α.

###### References

[Junker2008Analysis](@citet*)
"""
function centrality(::Type{KatzCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Union{Binary, Probabilistic}}; α::AbstractFloat=0.1)    
    @assert 0.0 <= α

    A = Matrix(N.edges.edges)
    C = (inv(I - α*A')-I)*ones(richness(N))

	return Dict(zip(species(N), C ./ sum(C)))
end

@testitem "We cannot use Katz centrality with a negative attenuation factor" begin
    nodes = Unipartite([:A, :B, :C])
    edges = Binary(Bool[0 1 1; 1 0 1; 0 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test_throws AssertionError centrality(KatzCentrality, N; α = -0.5)
end

@testitem "Katz centralities sum to one" begin
    nodes = Unipartite([:A, :B, :C])
    edges = Binary(Bool[0 1 1; 1 0 1; 0 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test sum(values(centrality(KatzCentrality, N; α = 0.2))) == 1.0
end

"""
    centrality(::Type{EigenvectorCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions})

Eigencentrality, corrected so that the centralities in the network sum to one.

The eigenvector centrality is calculated using Von Mises iteration, starting
from a random vector.

where A is the adjacency matrix of the graph, b is the vector of centralities,
At each iteration, the vector is updated to be A×b/|A×b|, and |⋅| is the norm.

The number of iterations is determined by the variable
`SpeciesInteractionNetworks.CENTRALITY_MAXITER`, and the algorithm usually
converges rapidly.

###### References

[Landau1895Zur](@citet*)
"""
function centrality(::Type{EigenvectorCentrality}, N::SpeciesInteractionNetwork{<:Unipartite, <:Interactions})
    b = rand(richness(N))
    b ./= sum(b)

    for _ in 1:SpeciesInteractionNetworks.CENTRALITY_MAXITER
        b1 = Array(N) * b
        n = LinearAlgebra.norm(b1)
        if isnan(n)
            @warn "Early return of the eigencetrality function, results may be approximate"
            break
        else
            for i in eachindex(b1)
                b[i] = b1[i]/n
            end
        end
    end

    b ./= sum(b)

    return Dict(zip(species(N), b ./ sum(b)))

end

@testitem "Eigenvector centralities sum to one" begin
    nodes = Unipartite([:A, :B, :C])
    edges = Binary(Bool[0 1 1; 1 0 1; 0 1 0])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test sum(values(centrality(EigenvectorCentrality, N))) ≈ 1.0
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
=#
