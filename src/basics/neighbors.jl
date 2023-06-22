"""
    successors(N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions}, sp::T) where {T}

The successors of a species in a bipartite network is the list of all species it
establishes a non-zero interaction with. For probabilistic networks, this
includes all species with a non-zero probability of interaction.

If the species is at the bottom of the network, or if the specis has no
successors, this method will retun an empty list of species, specifically
`Set{T}()`.
"""
function successors(N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions}, sp::T) where {T}
    if !(sp in species(N))
        throw(ArgumentError("The species $(sp) is not in the network"))
    end
    if sp in N.nodes.bottom
        return Set{T}()
    end
    succ_idx = findall(!iszero, N[sp,:])
    if isempty(succ_idx)
        return Set{T}()
    end
    return Set{T}(N.nodes.bottom[succ_idx])
end

"""
    successors(N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions}, sp::T) where {T}

The successors of a species in a unipartite network is the list of all species
it establishes a non-zero interaction with. For probabilistic networks, this
includes all species with a non-zero probability of interaction.

If the specis has no successors, this method will retun an empty list of
species, specifically `Set{T}()`.
"""
function successors(N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions}, sp::T) where {T}
    if !(sp in species(N))
        throw(ArgumentError("The species $(sp) is not in the network"))
    end
    succ_idx = findall(!iszero, N[sp,:])
    if isempty(succ_idx)
        return Set{T}()
    end
    return Set{T}(N.nodes.margin[succ_idx])
end

@testitem "We cannot look for successors of a species not in a network" begin
    edges = Binary(rand(Bool, (4, 3)))
    nodes = Bipartite([:A, :B, :C, :D], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test_throws ArgumentError successors(N, :X)
end

@testitem "Bottom-level species have no successors" begin
    edges = Binary(rand(Bool, (4, 3)))
    nodes = Bipartite([:A, :B, :C, :D], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test isempty(successors(N, :a))
end

@testitem "We can correctly identify successors in bipartite networks" begin
    M = [true true true; true true false; false false true]
    edges = Binary(M)
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test successors(N, :A) == Set([:a, :b, :c])
    @test successors(N, :B) == Set([:a, :b])
    @test successors(N, :C) == Set([:c])
end

@testitem "We can correctly identify successors in unipartite networks" begin
    M = [true true true; true true false; false false true]
    edges = Binary(M)
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test successors(N, :A) == Set([:A, :B, :C])
    @test successors(N, :B) == Set([:A, :B])
    @test successors(N, :C) == Set([:C])
end

"""
    predecessors(N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions}, sp::T) where {T}

The predecessors of a species in a bipartite network is the list of all species
it receives a non-zero interaction from. For probabilistic networks, this
includes all species with a non-zero probability of interaction.

If the species is at the top of the network, or if the specis has no
predecessors, this method will retun an empty list of species, specifically
`Set{T}()`.
"""
function predecessors(N::SpeciesInteractionNetwork{Bipartite{T}, <:Interactions}, sp::T) where {T}
    if !(sp in species(N))
        throw(ArgumentError("The species $(sp) is not in the network"))
    end
    if sp in N.nodes.top
        return Set{T}()
    end
    succ_idx = findall(!iszero, N[:,sp])
    if isempty(succ_idx)
        return Set{T}()
    end
    return Set{T}(N.nodes.top[succ_idx])
end


"""
    predecessors(N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions}, sp::T) where {T}

The predecessors of a species in a unipartite network is the list of all species
it receives a non-zero interaction from. For probabilistic networks, this
includes all species with a non-zero probability of interaction.

If the specis has no predecessors, this method will retun an empty list of
species, specifically `Set{T}()`.
"""
function predecessors(N::SpeciesInteractionNetwork{Unipartite{T}, <:Interactions}, sp::T) where {T}
    if !(sp in species(N))
        throw(ArgumentError("The species $(sp) is not in the network"))
    end
    succ_idx = findall(!iszero, N[:,sp])
    if isempty(succ_idx)
        return Set{T}()
    end
    return Set{T}(N.nodes.margin[succ_idx])
end

@testitem "We cannot look for predecessors of a species not in a network" begin
    edges = Binary(rand(Bool, (4, 3)))
    nodes = Bipartite([:A, :B, :C, :D], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test_throws ArgumentError predecessors(N, :X)
end

@testitem "Top-level species have no predecessors" begin
    edges = Binary(rand(Bool, (4, 3)))
    nodes = Bipartite([:A, :B, :C, :D], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test isempty(predecessors(N, :A))
end

@testitem "We can correctly identify predecessors in bipartite networks" begin
    M = [true true true; true true false; false false true]
    edges = Binary(M)
    nodes = Bipartite([:A, :B, :C], [:a, :b, :c])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test predecessors(N, :a) == Set([:A, :B])
    @test predecessors(N, :b) == Set([:A, :B])
    @test predecessors(N, :c) == Set([:A, :C])
end

@testitem "We can correctly identify predecessors in unipartite networks" begin
    M = [true true true; true true false; false false true]
    edges = Binary(M)
    nodes = Unipartite([:A, :B, :C])
    N = SpeciesInteractionNetwork(nodes, edges)
    @test predecessors(N, :A) == Set([:A, :B])
    @test predecessors(N, :B) == Set([:A, :B])
    @test predecessors(N, :C) == Set([:A, :C])
end

"""
    neighbors(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}) where {T}

The neighbors of a species is the list of both its successors and prde
"""
function neighbors(N::SpeciesInteractionNetwork{<:Partiteness{T}, <:Interactions}, sp::T) where {T}
    return union(predecessors(N,sp), successors(N,sp))
end