"""
    linearfilter(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}; α::Vector{T}=ones(4)) where {T <: AbstractFloat}



###### References

[Stock2017Linear](@citet)
"""
function linearfilter(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}; α::Vector{T}=ones(4)) where {T <: AbstractFloat}
    @assert length(α) == 4
    @assert all(α .>= 0.0 )
    α = α ./ sum(α)

    edges = Probabilistic(
        α[1].*N.edges.edges
            .+ α[2].*mean(N.edges.edges, dims=2)
            .+ α[3].*mean(N.edges.edges, dims=1)
            .+ α[4].*mean(N.edges.edges)
        )
    
    return SpeciesInteractionNetwork(copy(N.nodes), edges)

end

"""
    nullmodel(::Type{Degree}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

TODO
"""
nullmodel(::Type{Degree}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 1.0, 1.0, 0.0])
nullmodel(::Type{Connectance}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 0.0, 0.0, 1.0])
nullmodel(::Type{Generality}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 1.0, 0.0, 0.0])
nullmodel(::Type{Vulnerability}, N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = linearfilter(N; α=[0.0, 0.0, 1.0, 0.0])