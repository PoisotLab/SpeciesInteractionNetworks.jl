module TestOverlap
  using Test
  using EcologicalNetworks

  # Generate some data

  A = BipartiteNetwork([true false; true true; false true], [:a, :b, :c], [:d, :e])
  @test length(AJS(A)) == 2
  @test length(EAJS(A)) == 2
  @test AJS(A)[1].second == 0.5
  @test EAJS(A)[1].second == 0.5
  @test length(overlap(A)) == 2
  @test overlap(A)[1].second == 1

  B = BipartiteNetwork([true true true; true true false; true false false], [:a, :b, :c], [:d, :e, :f])
  @test AJS(B) == EAJS(B)
  @test length(AJS(B)) == 3
  @test length(overlap(B)) == 3
  @test filter(x -> x.first == Set((:a, :b)), overlap(B))[1].second == 2

end
