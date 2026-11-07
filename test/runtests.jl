using Test
using BioExplorer
using DataFrames
using Random

@testset "Darwin Core ingestion" begin
    records = DataFrame(
        eventID = ["A", "A", "B"],
        scientificName = ["x", "y", "x"],
        individualCount = [1, 2, 3],
    )

    community, diagnostics, exclusions = from_records(records)

    @test community.values == [1.0 1.0; 1.0 0.0]
    @test isempty(exclusions)

    @test_throws ArgumentError from_records(
        records;
        semantics = :abundance,
    )

    abundance, _, abundance_exclusions = from_records(
        records;
        semantics = :abundance,
        count_col = :individualCount,
    )

    @test abundance.values == [1.0 2.0; 3.0 0.0]
    @test isempty(abundance_exclusions)

    @test_throws ArgumentError CommunityMatrix(
        ["a"],
        ["x"],
        reshape([-1], 1, 1),
    )
end

@testset "Analyses" begin
    community = CommunityMatrix(
        ["A", "B"],
        ["x", "y"],
        [1.0 1.0; 1.0 0.0];
        semantics = :incidence,
    )

    alpha = alpha_diversity(community)
    beta = beta_diversity(community)

    @test alpha.shannon[1] ≈ log(2)
    @test beta.jaccard[1] ≈ 0.5
    @test beta.sorensen[1] ≈ 1 / 3

    original_values = copy(community.values)

    curve1 = accumulation_curve(
        community;
        rng = MersenneTwister(1),
    )
    curve2 = accumulation_curve(
        community;
        rng = MersenneTwister(1),
    )

    @test curve1 == curve2
    @test community.values == original_values

    # Incidence data should produce Chao2 output.
    estimators = richness_estimators(community)

    @test names(estimators) == [
        "observed_richness",
        "uniques",
        "duplicates",
        "sites_sampled",
        "chao2",
    ]
    @test nrow(estimators) == 1

    # Abundance data should produce Chao1 output.
    abundance = CommunityMatrix(
        ["a"],
        ["x", "y"],
        [1 2];
        semantics = :abundance,
    )

    abundance_estimators = richness_estimators(abundance)

    @test "chao1" in names(abundance_estimators)
    @test abundance_estimators.chao1[1] ≈ 2.5
    @test sampling_coverage(abundance).sample_coverage[1] ≈ 2 / 3

    empty_community = CommunityMatrix(
        ["empty"],
        ["x"],
        zeros(1, 1);
        semantics = :incidence,
    )

    @test alpha_diversity(empty_community).hill_q1[1] == 0
end