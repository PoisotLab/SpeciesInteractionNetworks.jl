"""
    null4(N::BinaryNetwork)

Given a matrix `A`, `null4(A)` returns a matrix with the same dimensions, where
every interaction happens with a probability equal to the product of the degree
of each species.
"""
function null4(N::BinaryNetwork)
  Afiltered = sum(N, dims=1) .* sum(N, dims=2) ./ sum(N)^2
  ReturnType = typeof(N) <: AbstractBipartiteNetwork ? BipartiteProbabilisticNetwork : UnipartiteProbabilisticNetwork
  return ReturnType(Afiltered, _species_objects(N)...)
end
