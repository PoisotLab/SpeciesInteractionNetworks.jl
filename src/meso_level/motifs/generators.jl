"""
    motifs(::Type{Bipartite}, nodes::Integer)

Returns a list of bipartite motifs with the specified number of nodes. The list
of motifs for each richness are from [Simmons2018Motifs](@citet).

###### References

[Simmons2018Motifs](@citet*)
"""
function motifs(::Type{Bipartite}, nodes::Integer)
    A = Matrix{Bool}[]
    if nodes == 2
        A = [ones(Bool, 1,1)]
    end
    if nodes == 3
        A = Matrix{Bool}[
            Bool[1 1;],
            reshape(Bool[1; 1], (2,1))
        ]
    end
    if nodes == 4
        A = Matrix{Bool}[
            Bool[1 1 1;],
            Bool[1 1; 0 1],
            Bool[1 1; 1 1],
            reshape(Bool[1; 1; 1], (3,1)),
        ]
    end
    if nodes == 5
        A = Matrix{Bool}[
            Bool[1 1 1 1;],
            Bool[1 0; 1 0; 1 1;],
            Bool[1 0; 1 1; 0 1;],
            Bool[1 0; 1 1; 1 1;],
            Bool[1 1; 1 1; 1 1;],
            Bool[1 0 0; 1 1 1;],
            Bool[1 1 0; 0 1 1;],
            Bool[1 1 0; 1 1 1;],
            Bool[1 1 1; 1 1 1],
            reshape(Bool[1; 1; 1; 1;], (4, 1)),
        ]
    end
    if isempty(A)
        throw("Bipartite motifs with $(nodes) nodes are not supported")
    end
    return SpeciesInteractionNetwork{Bipartite, Binary}.(A)
end

"""
    motifs(::Type{Unipartite}, nodes::Integer)

Returns a list of unipartite motifs with the specified number of nodes. The list
of motifs for three species are from [Stouffer2007Evidence](@citet).

###### References

[Stouffer2007Evidence](@citet*)
"""
function motifs(::Type{Unipartite}, nodes::Integer)
    A = Matrix{Bool}[]
    if nodes == 2
        A = [ones(Bool, 1,1)]
    end
    if nodes == 3
        A = Matrix{Bool}[
            Bool[0 1 0; 0 0 1; 0 0 0],
            Bool[0 1 1; 0 0 1; 0 0 0],
            Bool[0 1 0; 0 0 1; 1 0 0],
            Bool[0 0 1; 0 0 1; 0 0 0],
            Bool[0 1 1; 0 0 0; 0 0 0],
            Bool[0 1 1; 0 0 0; 1 1 0],
            Bool[0 1 1; 0 0 1; 0 1 0],
            Bool[0 0 1; 0 0 0; 1 1 0],
            Bool[0 1 0; 0 0 1; 0 1 0],
            Bool[0 1 0; 0 0 1; 1 1 0],
            Bool[0 1 1; 1 0 1; 1 1 0],
            Bool[0 1 1; 1 0 0; 1 1 0],
            Bool[0 1 1; 1 0 0; 1 0 0]
        ]
    end
    if isempty(A)
        throw("Unipartite motifs with $(nodes) nodes are not supported")
    end
    return SpeciesInteractionNetwork{Unipartite, Binary}.(A)
end