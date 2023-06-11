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
    
function swap!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}, constraint::Type{PC}) where {PC <: PermutationConstraint}
    if length(N) < 2 
        throw(ArgumentError("Impossible to swap a network with fewer than two interactions"))
    end
    (constraint == Degree) && swap_degree!(N)
    return N
end

@testitem "We cannot swap a network with fewer than two interactions" begin
    nodes = Unipartite([:a, :b, :c])
    edges = Binary(zeros(Bool, (3,3)))
    N = SpeciesInteractionNetwork(nodes, edges)
    @test_throws ArgumentError swap!(N, Degree)
end

function swap_degree!(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})
    intpool = interactions(N)
    intpair = StatsBase.sample(intpool, 2, replace=false)
    swappair = [(intpair[1][1],intpair[2][2],true), (intpair[2][1],intpair[1][2],true)]
    valid = all(iszero.(N, swappair))
    iters = 0
    while (!valid)||(iter < SWAP_MAXITER)
        intpair = StatsBase.sample(intpool, 2, replace=false)
        swappair = [(intpair[1][1],intpair[2][2],true), (intpair[2][1],intpair[1][2],true)]
        valid = all(iszero.(N, swappair))
        iters += 1
    end
    if iter < SWAP_MAXITER
        N[intpair[1]] = false
        N[intpair[2]] = false
        N[swappair[1]] = true
        N[swappair[2]] = true
    end
    return N
end