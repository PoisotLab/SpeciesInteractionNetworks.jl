LinearAlgebra.svd(N::SpeciesInteractionNetwork) = LinearAlgebra.svd(Array(N))
LinearAlgebra.rank(N::SpeciesInteractionNetwork) = LinearAlgebra.rank(Array(N))

"""
    tsvd(N::SpeciesInteractionNetwork, r::Integer)

Returns the left and right subspaces for a given ecological network, applying a
truncation of the SVD at rank `r`. Note that if `r` is higher than the rank of
the network, it will be scaled appropriately.

###### References

[DallaRiva2016Exploring](@citet*)
"""
function tsvd(N::SpeciesInteractionNetwork, r::Integer)
    r = min(r, rank(N))
    U, Σ, V = svd(N)
    L = similar(U)
    R = similar(V)
    mul!(L, U, diagm(sqrt.(Σ)))
    mul!(R, diagm(sqrt.(Σ)), V')
    L = L[:,1:r]
    R = R[1:r,:]
    return L, R
end

"""
    rdpg(N::SpeciesInteractionNetwork, r::Integer)

Returns a quantitative network predicted using a random dot-product graph on the
[`tsvd`](@ref) of the network at rank `r`.

RDPG has been shown to be predictive of food web structure even in the context
of transfer learning [Strydom2022Food](@cite).

###### References

[DallaRiva2016Exploring](@citet*)

[Strydom2022Food](@citet*)
"""
function rdpg(N::SpeciesInteractionNetwork, r::Integer)
    L, R = tsvd(N, r)
    A = clamp.(L * R, zero(eltype(L)), one(eltype(L)))
    return SpeciesInteractionNetwork(copy(N.nodes), Quantitative(A))
end

"""
    rdpg(N::SpeciesInteractionNetwork, r::Integer, threshold::T) where {T <: AbstractFloat}

Returns a binary network predicted using a random dot-product graph on the
[`tsvd`](@ref) of the network at rank `r`, where interactions with a score above
`threshold` are true.

This approach to RDPG has been shown to be highly predictive under strong data sparsity [Poisot2023Network](@cite).

###### References

[DallaRiva2016Exploring](@citet*)

[Poisot2023Network](@citet*)
"""
function rdpg(N::SpeciesInteractionNetwork, r::Integer, threshold::T) where {T <: AbstractFloat}
    L, R = tsvd(N, r)
    A = (L*R).>= threshold
    return SpeciesInteractionNetwork(copy(N.nodes), Binary(A))
end

"""
    complexity(N::SpeciesInteractionNetwork)

Returns the SVD complexity of a network, defined as the Pielou entropy of its
singular values. Note that only the singular values up to the rank of the
network are returned.

###### References

[Strydom2021SVD](@citet*)
"""
function complexity(N::SpeciesInteractionNetwork)
    U, Σ, V = svd(N)
    Σ = Σ[1:rank(N)]
    p = Σ ./ sum(Σ)
    H = -sum(p .* log.(p))
    return H / log(rank(N))
end