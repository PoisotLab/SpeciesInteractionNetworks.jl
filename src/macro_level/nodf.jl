"""
WNODF of a single axis
"""
function nodf_axis(N::SpeciesInteractionNetwork{<:Bipartite, <:Quantitative})
    A = Array(N)
    row_order = sortperm(vec(sum(A; dims = 2)); rev = true)
    A = map(Float64, A[row_order, :])

    WNODFr = 0.0
    for i in 1:(size(A, 1) - 1)
        for j in (i + 1):size(A, 1)
            if (sum(A[i, :]) .> sum(A[j, :]))
                Nj = sum(A[j, :] .> 0.0)
                kij = 0.0
                # The following is only applied to interactions that are non-0
                # in j, which is not really clear in the original paper
                for l in eachindex(A[i, :])
                    if (A[j, l] > 0.0)
                        if A[j, l] < A[i, l]
                            kij += 1.0
                        end
                    end
                end
                WNODFr += kij / Nj
            end
        end
    end

    correction = (richness(N, 1) * (richness(N, 1) - 1))
    return 2WNODFr / correction
end

"""
NODF of a single axis
"""
function nodf_axis(N::SpeciesInteractionNetwork{<:Bipartite, <:Binary})
    A = Array(N)
    row_order = sortperm(vec(sum(A; dims = 2)); rev = true)
    A = A[row_order, :]
    d = float(vec(sum(A; dims = 2)))

    NODFr = 0.0
    for i in 1:(size(A, 1) - 1)
        for j in (i + 1):size(A, 1)
            DFpaired = d[j] < d[i] ? 1.0 : 0.0
            Npaired = sum(A[i, :] .& A[j, :]) / d[j]
            NODFr += (DFpaired * Npaired)
        end
    end
    correction = (richness(N, 1) * (richness(N, 1) - 1))
    return 2NODFr / correction
end

"""
    nodf(N::T; dims::Union{Nothing,Integer}=nothing) where {T <: Union{BipartiteNetwork,BipartiteQuantitativeNetwork}}

Returns `nodf` for a margin of the network. The `i` argument can be 1 for
top-level, 2 for bottom level, and the function will throw an `ArgumentError` if
an invalid value is used. For quantitative networks, *WNODF* is used.

###### References

[AlmeidaNeto2008consistent](@citet*)

[AlmeidaNeto2011straightforward](@citet*)
"""
function nodf(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Quantitative}})
    return (nodf(N, 1) + nodf(N, 2)) / 2
end

function nodf(N::SpeciesInteractionNetwork{<:Bipartite, <:Union{Binary, Quantitative}}, dims::Integer)
    if dims == 1
        return nodf_axis(N)
    end
    if dims == 2
        return nodf_axis(permutedims(N))
    end
    throw(
        ArgumentError(
            "dims can only be 1 (nestedness of rows) or 2 (nestedness of columns), you used $(dims)",
        ),
    )
end
