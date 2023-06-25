"""
    labelpropagation(N::SpeciesInteractionNetwork)

Uses label propagation to generate a first approximation of the modular
structure of a network. This is usually followed by the BRIM (`brim`) method.
This method supposedly performs better for large graphs, but we rarely observed
any differences between it and variations of BRIM alone on smaller graphs.

###### References

[Liu2010Community](@citet*)
"""
function labelpropagation(N::SpeciesInteractionNetwork)
    L = collect(1:richness(N))
    
    m₀ = modularity(BRIM, N, L)
    mₜ = m₀
    
    improvement = true
    cursor = 0

    while (improvement)&(cursor < MODULARITY_MAXITER)
        cursor += 1
        order_1 = Random.shuffle(1:richness(N,1))
        order_2 = Random.shuffle(1:richness(N,2))

        mₜ = modularity(BRIM, N, L)
        m₀, improvement = mₜ > m₀ ? (mₜ, true) : (mₜ, false)
    end

end