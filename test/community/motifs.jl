module TestMotifs
using Test
#include("../../src/SpeciesInteractionNetworks.jl")
using EcologicalNetworks
using LinearAlgebra

@test unipartitemotifs().S1[:a, :b]
@test !unipartitemotifs().S1[:a, :a]

n = UnipartiteNetwork(zeros(Bool, (2, 2)))
m = unipartitemotifs().S1
@test length(find_motif(n, m)) ≈ 0.0

# Test with a single link
N = UnipartiteProbabilisticNetwork([0.2 0.8; 0.2 0.1])
m = UnipartiteNetwork([false true; false false])
@test expected_motif_count(find_motif(N,m))[1] ≈ 0.680

# Test with a perfect match
N = UnipartiteProbabilisticNetwork([0.0 1.0; 0.0 0.0])
m = UnipartiteNetwork([false true; false false])
@test expected_motif_count(find_motif(N,m)) == (1.0, 0.0)

# Test with a fork-like thing
diam = UnipartiteNetwork([0 1 0 0; 0 0 1 1; 0 0 0 0; 0 0 0 0].>0)
clin = UnipartiteNetwork([0 1 0; 0 0 1; 0 0 0].>0)
capp = UnipartiteNetwork([0 1 1; 0 0 0; 0 0 0].>0)
cdir = UnipartiteNetwork([0 0 1; 0 0 1; 0 0 0].>0)

@test length(find_motif(diam, clin)) == 2
@test length(find_motif(diam, capp)) == 1
@test length(find_motif(diam, cdir)) == 0

# Test with a diamond food web
diam = UnipartiteNetwork([0 1 1 0; 0 0 0 1; 0 0 0 1; 0 0 0 0].>0)
clin = UnipartiteNetwork([0 1 0; 0 0 1; 0 0 0].>0)
capp = UnipartiteNetwork([0 1 1; 0 0 0; 0 0 0].>0)
cdir = UnipartiteNetwork([0 0 1; 0 0 1; 0 0 0].>0)

@test length(find_motif(diam, clin)) == 2
@test length(find_motif(diam, capp)) == 1
@test length(find_motif(diam, cdir)) == 1

lfc = UnipartiteNetwork([0 1 0; 0 0 1; 0 0 0].>0)
wsl = UnipartiteNetwork([1 1 0; 0 0 1; 0 0 0].>0)

# This tests for the removal of the diagonal
@test length(find_motif(wsl, lfc)) == 1

lfc = UnipartiteNetwork([0 1 0; 0 0 1; 0 0 0].>0)
wsl = UnipartiteQuantitativeNetwork([1 0.5 0; 0 0 1.6; 0 0 0])

@test length(find_motif(wsl, lfc)) == 1

# Test with a diamond food web
diam = BipartiteNetwork([1 1 0; 1 1 1].>0)
tthr = BipartiteNetwork([1 0; 1 1].>0)
tfou = BipartiteNetwork([1 1; 1 1].>0)

@test length(find_motif(diam, tthr)) == 2
@test length(find_motif(diam, tfou)) == 1

# Test on a three-species network
B = UnipartiteNetwork([false true true; false false true; false false false])
@test length(find_motif(B, B)) == 1

BDN = BipartiteNetwork([false true true; false false true; false false false])
@test length(find_motif(BDN, BDN)) == 1.0

# Test on the same network, this time with a proba one
P = UnipartiteProbabilisticNetwork([0.0 1.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0])
B = UnipartiteNetwork([false true true; false false true; false false false])
@test expected_motif_count(find_motif(P, B)) == (1.0, 0.0)

# Test with some variation
ovl = UnipartiteNetwork([false true true; false false true; false false false])
N = UnipartiteProbabilisticNetwork([0.0 1.0 0.8; 0.0 0.0 1.0; 0.0 0.0 0.0])
@test first(expected_motif_count(find_motif(N, ovl))) ≈ 0.8
fchain = UnipartiteNetwork([false true false; false false true; false false false])
@test first(expected_motif_count(find_motif(N, fchain))) ≈ 0.2

# Bipartite proba network
B = BipartiteProbabilisticNetwork(Matrix{Float64}(I, (3, 3)))
m = BipartiteNetwork(Matrix{Bool}(I, (3, 3)))
@test first(expected_motif_count(find_motif(B, m))) ≈ 1.0
@test last(expected_motif_count(find_motif(B, m))) ≈ 0.0

# Unipartite multiple motifs
X = [false true true true false; false false false true true; false false false false true;
    false false false false false; false false false false false]
N = UnipartiteNetwork(X, [:a, :b, :c, :d, :e])

@test length(find_motif(N, unipartitemotifs().S1)) == 2
@test length(find_motif(N, unipartitemotifs().S2)) == 1
@test length(find_motif(N, unipartitemotifs().S3)) == 0
@test length(find_motif(N, unipartitemotifs().S4)) == 1
@test length(find_motif(N, unipartitemotifs().S5)) == 3

# Diamond
N = UnipartiteNetwork([false true true false; false false false true;
    false false false true; false false false false], [:a, :b, :c, :d])

@test length(find_motif(N, unipartitemotifs().S1)) == 2
@test length(find_motif(N, unipartitemotifs().S2)) == 0
@test length(find_motif(N, unipartitemotifs().S3)) == 0
@test length(find_motif(N, unipartitemotifs().S4)) == 1
@test length(find_motif(N, unipartitemotifs().S5)) == 1

# Self motifs
for this_motif in unipartitemotifs()
    @test length(find_motif(this_motif, this_motif)) == 1
end

end
