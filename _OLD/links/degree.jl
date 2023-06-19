"""
**Variance in the outgoing degree**

    degree_out_var(N::ProbabilisticNetwork)
"""
function degree_out_var(N::ProbabilisticNetwork)
  d_o_v = vec(mapslices(_additive_bernoulli_variance, convert(Matrix, N.edges); dims=2))
  return Dict(zip(species(N; dims=1), d_o_v))
end

"""
**Variance in the ingoing degree**

    degree_in_var(N::ProbabilisticNetwork)
"""
function degree_in_var(N::ProbabilisticNetwork)
  d_i_v = vec(mapslices(_additive_bernoulli_variance, convert(Matrix, N.edges); dims=1))
  return Dict(zip(species(N; dims=2), d_i_v))
end

"""
    degree_var(N::ProbabilisticNetwork; dims::Union{Nothing,Integer}=nothing)

Variance in the degree of species in a probabilistic network.
"""
function degree_var(N::ProbabilisticNetwork; dims::Union{Nothing,Integer}=nothing)
  dims == 1 && return degree_out_var(N)
  dims == 2 && return degree_in_var(N)
  if isnothing(dims)
    if typeof(N) <: AbstractBipartiteNetwork
      return merge(degree_out_var(N), degree_in_var(N))
    else
      din = degree_in_var(N)
      dout = degree_out_var(N)
      for k in keys(din)
        din[k] += dout[k]
      end
      return din
    end
  end
  throw(ArgumentError("dims can only be 1 (out-degre) or 2 (in-degree) or `nothing` (both), you used $(dims)"))
end
