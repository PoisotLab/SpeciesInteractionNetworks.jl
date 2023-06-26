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
    
    partition = Dict([species(N)[i]=>i for i in 1:richness(N)])
    
    m₀ = modularity(N, partition)
    mₜ = m₀
    
    improvement = true
    cursor = 0

    while (improvement)&(cursor < MODULARITY_MAXITER)
        cursor += 1
        order_1 = Random.shuffle(species(N, 1))
        order_2 = Random.shuffle(species(N, 2))

        for s1 in order_1
            succ = collect(successors(N, s1))
            labels = [partition[s] for s in succ]
            wts = [Float32(N[s1,s]) for s in succ]
            new_label = last(findmax(StatsBase.countmap(labels, wts)))
            partition[s1] = new_label
        end

        for s2 in order_2
            pred = collect(predecessors(N, s2))
            labels = [partition[s] for s in pred]
            wts = [Float32(N[s,s2]) for s in pred]
            new_label = last(findmax(StatsBase.countmap(labels, wts)))
            partition[s2] = new_label
        end

        mₜ = modularity(N, partition)
        m₀, improvement = mₜ > m₀ ? (mₜ, true) : (mₜ, false)
    end

    return partition

end