LinearAlgebra.svd(N::SpeciesInteractionNetwork) = LinearAlgebra.svd(Array(N))
LinearAlgebra.rank(N::SpeciesInteractionNetwork) = LinearAlgebra.rank(Array(N))

"""
    rdpg(N::SpeciesInteractionNetwork, r::Integer=3)

Returns the random dot-product graph left and right subspaces for a given
ecological network, applying a truncation at rank `r`. Note that if `r` is
higher than the rank of the network, it will be scaled appropriately.

RDPG has been shown to be predictive of food web structure even in the context
of transfer learning [Strydom2022Food](@cite) or extreme data sparsity
[Poisot2023Network](@cite).

###### References

[DallaRiva2016Exploring](@citet*)

[Poisot2023Network](@citet*)

[Strydom2022Food](@citet*)
"""
function rdpg(N::SpeciesInteractionNetwork, r::Integer=3)
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