using Test, BioExplorer, DataFrames, Random
@testset "Darwin Core ingestion" begin
 d=DataFrame(eventID=["A","A","B"],scientificName=["x","y","x"],individualCount=[1,2,3])
 c,diag,exc=from_records(d)
 @test c.values==[1 1;1 0]
 @test isempty(exc)
 @test_throws ArgumentError from_records(d;semantics=:abundance)
 a,_,_=from_records(d;semantics=:abundance,count_col=:individualCount)
 @test a.values==[1 2;3 0]
 @test_throws ArgumentError CommunityMatrix(["a"],["x"],reshape([-1],1,1))
end
@testset "Analyses" begin
 c=CommunityMatrix(["A","B"],["x","y"],[1.0 1;1 0];semantics=:incidence)
 @test alpha_diversity(c).shannon[1]≈log(2)
 @test beta_diversity(c).jaccard[1]≈0.5
 @test beta_diversity(c).sorensen[1]≈1/3
 x=copy(c.values)
 @test accumulation_curve(c;rng=MersenneTwister(1))==accumulation_curve(c;rng=MersenneTwister(1))
 @test x==c.values
 @test_throws ArgumentError richness_estimators(c)
 a=CommunityMatrix(["a"],["x","y"],[1 2])
 @test richness_estimators(a).chao1[1]==2.5
 @test sampling_coverage(a).sample_coverage[1]≈2/3
 e=CommunityMatrix(["empty"],["x"],zeros(1,1))
 @test alpha_diversity(e).hill_q1[1]==0
end
