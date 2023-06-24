using Revise
using SpeciesInteractionNetworks
import DelimitedFiles
using Statistics

data_path = joinpath(@__DIR__, "data", "nz_stream")
netfiles = filter(endswith(".csv"), readdir(data_path; join=true))

function readnz(nf)
    X = DelimitedFiles.readdlm(nf, ',')
    A = permutedims(X[2:end,2:end])
    S = [replace(s, r"^(.+) \(N=\d+\)$" => s"\1") for s in X[1,2:end]]
    nodes = Unipartite(S)
    edges = Binary(Bool.(A))
    return SpeciesInteractionNetwork(nodes, edges)
end

N = readnz(netfiles[1])
M = copy(N)

c0 = complexity(N)

Δ = c0 - complexity(M)

T0 = 0.01
λ = 0.99
budget = 5000

for t in 1:budget
    @info t
    Tt = T0*(λ^t)
    # Generate a candidate neighbor
    i_before = interactions(M)
    swap!(M, Degree)
    i_after = interactions(M)
    # Get the complexity of the candidate
    ct = complexity(M)
    # Get the difference
    Δ = c0 - ct
    # Decide
    if Δ < 0
        c0 = ct
        @info " → valid move"
    else
        P = exp(-(Δ)/Tt)
        if rand() <= P
            c0 = ct
            @info " → allowed move"
        else
            @info " → rejected move"
            intbefore = setdiff(i_before, i_after)
            intafter = setdiff(i_after, i_before)
            if !isempty(intbefore)
                for ib in intbefore
                    M[ib[1],ib[2]] = ib[3]
                end
                for ia in intafter
                    M[ia[1],ia[2]] = !ia[3]
                end
            end
        end
    end
    @info " → $(c0)"
end