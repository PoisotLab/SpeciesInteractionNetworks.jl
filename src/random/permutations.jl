const SWAP_MAXITER = 100

"""
    PermutationConstraint

The `PermutationConstraint` specifies which structural constraint is enforced.
It is defined as an abstract type, and the subtypes can be passed as the second
argument to `swap!`.

Currently supported constraints are `Degree` (degree distribution is
maintained), `Generality` (number of out-going links is maintained),
`Vulnerability` (number of in-going links is maintained), and `Connectance`
(only the connectance is maintained). Note that *in addition*, species cannot
become disconnected, even if the constraint is not acting on the degree / degree
distribution.
"""
abstract type PermutationConstraint end
abstract type Degree <: PermutationConstraint end
abstract type Generality <: PermutationConstraint end
abstract type Vulnerability <: PermutationConstraint end
abstract type Connectance <: PermutationConstraint end

"""
    swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

Performs *one* swap of interactions in the network. If no
`PermutationConstraint` is given as a second argument, the degree distribution
of all species will be maintained.
"""
swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}) = swap!(N, Degree)

"""
    swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Degree})

Permutations with a constraint by degree work by picking two interacting species
pairs, (r1, s1) and (r2, s2), and trying to replace them by (r1, s2) and (r2,
s1).
"""
function swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Degree})
    if length(N) < 2 
        throw(ArgumentError("Impossible to swap a network with fewer than two interactions"))
    end
    intpool = interactions(N)
    intpair = StatsBase.sample(intpool, 2, replace=false)
    swappair = [(intpair[1][1],intpair[2][2],true), (intpair[2][1],intpair[1][2],true)]
    valid = all(iszero.([N[swp] for swp in swappair]))
    iters = 0
    while (!valid)&&(iters < SpeciesInteractionNetworks.SWAP_MAXITER)
        intpair = StatsBase.sample(intpool, 2, replace=false)
        swappair = [(intpair[1][1],intpair[2][2],true), (intpair[2][1],intpair[1][2],true)]
        valid = all(iszero.([N[swp] for swp in swappair]))
        iters += 1
    end
    if iters <= SpeciesInteractionNetworks.SWAP_MAXITER
        N[intpair[1]] = false
        N[intpair[2]] = false
        N[swappair[1]] = true
        N[swappair[2]] = true
    end
    return N
end

"""
    swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Generality})

Permutations with a constraint by degree work by picking one interacting species
pair, (r1, s1), and a new stem species s3, trying to replace them by (r1, s3).
This function only applies if the result of this permutations does not remove
the last incoming link from s1.
"""
function swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Generality})
    if length(N) < 2 
        throw(ArgumentError("Impossible to swap a network with fewer than two interactions"))
    end
    intpool = interactions(N)
    intpair = rand(intpool)
    newstem = rand(species(N, 2))
    valid = (length(predecessors(N, intpair[2])) > 1)&(iszero(N[intpair[1], newstem]))
    iters = 0
    while (!valid)&&(iters < SpeciesInteractionNetworks.SWAP_MAXITER)
        intpair = rand(intpool)
        newstem = rand(species(N, 2))
        valid = (length(predecessors(N, intpair[2])) > 1)&(iszero(N[intpair[1], newstem]))
        iters += 1
    end
    if iters <= SpeciesInteractionNetworks.SWAP_MAXITER
        N[intpair] = false
        N[intpair[1], newstem] = true
    end
    return N
end

"""
    swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Vulnerability})

Permutations with a constraint by degree work by picking one interacting species
pair, (r1, s1), and a new root species r3, trying to replace them by (r3, s1).
This function only applies if the result of this permutations does not remove
the last outgoing link from r1.
"""

function swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Vulnerability})
    if length(N) < 2 
        throw(ArgumentError("Impossible to swap a network with fewer than two interactions"))
    end
    intpool = interactions(N)
    intpair = rand(intpool)
    newroot = rand(species(N, 1))
    valid = (length(successors(N, intpair[1])) > 1)&(iszero(N[newroot, intpair[2]]))
    iters = 0
    while (!valid)&&(iters < SpeciesInteractionNetworks.SWAP_MAXITER)
        intpair = rand(intpool)
        newroot = rand(species(N, 1))
        valid = (length(successors(N, intpair[1])) > 1)&(iszero(N[newroot, intpair[2]]))
        iters += 1
    end
    if iters <= SpeciesInteractionNetworks.SWAP_MAXITER
        N[intpair] = false
        N[newroot, intpair[2]] = true
    end
    return N
end

"""
    swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Connectance})

Permutations with a constraint by connectance will *randomly* (and with equal
probability) perform a move that is constrained by degree, generality, or
vulnerability.
"""
function swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, ::Type{Connectance})
    strategy = rand([Generality, Vulnerability, Degree])
    swap!(N, strategy)
    return N
end

@testitem "We can degree-swap a unipartite network with enough interactions" begin
    nodes = Unipartite([:A, :B, :C, :D, :E, :F])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    it_orig = interactions(N)
    for i in 1:10
        swap!(N, Degree)
    end
    it_swap = interactions(N)
    @test it_orig !== it_swap
end

@testitem "We can degree-swap a bipartite network with enough interactions" begin
    nodes = Bipartite([:A, :B, :C, :D, :E, :F], [:a, :b, :c, :d, :e, :f, :g, :h])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    it_orig = interactions(N)
    for i in 1:10
        swap!(N, Degree)
    end
    it_swap = interactions(N)
    @test it_orig !== it_swap
end

@testitem "We can generality-swap a bipartite network with enough interactions" begin
    nodes = Bipartite([:A, :B, :C, :D, :E, :F], [:a, :b, :c, :d, :e, :f, :g, :h])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    it_orig = interactions(N)
    for i in 1:10
        swap!(N, Generality)
    end
    it_swap = interactions(N)
    @test it_orig !== it_swap
end

@testitem "We can vulnerability-swap a bipartite network with enough interactions" begin
    nodes = Bipartite([:A, :B, :C, :D, :E, :F], [:a, :b, :c, :d, :e, :f, :g, :h])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    it_orig = interactions(N)
    for i in 1:10
        swap!(N, Vulnerability)
    end
    it_swap = interactions(N)
    @test it_orig !== it_swap
end

@testitem "We can connectance-swap a bipartite network with enough interactions" begin
    nodes = Bipartite([:A, :B, :C, :D, :E, :F], [:a, :b, :c, :d, :e, :f, :g, :h])
    edges = Binary(rand(Bool, (richness(nodes,1), richness(nodes, 2))))
    N = SpeciesInteractionNetwork(nodes, edges)
    it_orig = interactions(N)
    for i in 1:10
        swap!(N, Connectance)
    end
    it_swap = interactions(N)
    @test it_orig !== it_swap
end

@testitem "We cannot swap a network with fewer than two interactions" begin
    nodes = Unipartite([:a, :b, :c])
    edges = Binary(zeros(Bool, (3,3)))
    N = SpeciesInteractionNetwork(nodes, edges)
    @test_throws ArgumentError swap!(N, Degree)
end
