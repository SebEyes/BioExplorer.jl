"""
    alpha_diversity(c::CommunityMatrix)

Calculate basic alpha-diversity metrics for each site in a community matrix.

The function processes each row of `c.values` independently and returns a
`DataFrame` with one row per site. Taxa with values greater than zero are
treated as observed. The function does not modify `c`.

The input matrix should contain non-negative values. For abundance-based
interpretation, values should represent comparable taxon abundances or counts.
When the matrix uses incidence data, the resulting proportional metrics treat
each observed taxon equally.

# Argument
- `c::CommunityMatrix`: Validated site-by-taxon matrix. Rows correspond to
  sites, columns correspond to taxa, and `c.sites` supplies the site
  identifiers used in the output.

# Returned columns
- `site`: Site identifier.
- `richness`: Number of taxa with a value greater than zero.
- `shannon`: Shannon entropy, calculated as
  `-sum(pᵢ * log(pᵢ))` over taxa with positive proportions.
- `simpson`: Simpson concentration, calculated as `sum(pᵢ^2)`. This is the
  concentration form of the Simpson index, not `1 - sum(pᵢ^2)`.
- `hill_q1`: Order-1 Hill number, calculated as `exp(shannon)`.
- `hill_q2`: Order-2 Hill number, calculated as `1 / simpson`.
- `pielou`: Pielou's evenness, calculated as `shannon / log(richness)`.
  Returns `missing` when richness is zero or one, for which this expression
  is undefined.

Here, `pᵢ` is the value for taxon `i` divided by the sum of values at that site.
If a site's total is zero, its proportions are treated as zero; its Shannon,
Simpson, and Hill-number results are consequently zero.

# Notes
The metrics describe the values recorded in the matrix. A zero in the matrix
means no value was recorded for that site–taxon pair; it does not, by itself,
establish a verified ecological absence.

# Example
```julia
results = alpha_diversity(community)
```
"""
function alpha_diversity(c::CommunityMatrix)
    # Initialise output dataframe
    d=DataFrame(
        site=String[],                      # Site ID
        richness=Int[],                     # Taxonomic richness (i.e. number of recorded taxon)
        shannon=Float64[],                  # Shannon index
        simpson=Float64[],                  # Simpson index
        hill_q1=Float64[],                  # Hill number Q = 1
        hill_q2=Float64[],                  # Hill number Q = 2
        pielou=Union{Missing,Float64}[]     # Pielou eveness
    )

    # Comput alpha diversity metrics for all sites recursively
    for (i,s) in enumerate(c.sites)
        x=Float64.(c.values[i,:]);                  # Collect values (incidence or abundance)
        n=sum(x);                                   # Total collected abundance
        p=n==0 ? zeros(length(x)) : x./n            # Proportion of each taxon
        r=count(>(0),x);                            # Taxonomic richness
        h=-sum(q>0 ? q*log(q) : 0.0 for q in p);    # Shannon index
        si=sum(p.^2)                                # Simpson index

        push!(
            d,
            (
                s,                                  # Site ID
                r,                                  # Taxonomic richness
                h,
                si,
                r==0 ? 0.0 : exp(h),                # Hill q=1 (Exponential Shannon)
                si==0 ? 0.0 : 1/si,                 # Hill q = 2 (Inverse simpson)
                r<=1 ? missing : h/log(r)           # Pielou eveness
            )
        )
    end

    # Return output as DataFrame
    d
end

"""
    beta_diversity(c::CommunityMatrix)

Compute pairwise beta-diversity dissimilarities among all sites in a
`CommunityMatrix`.

Each unordered pair of sites is compared once. The function returns a
`DataFrame` with one row per pair and columns for Jaccard, Sørensen and
Bray–Curtis dissimilarities, together with total dissimilarity and its
turnover and richness-difference components according to Baselga and
Carvalho.

# Argument
- `c::CommunityMatrix`: A validated site-by-taxon matrix. Rows represent
  sites and columns represent taxa. Values must be finite and non-negative,
  as enforced by the `CommunityMatrix` constructor.

# Returned value
A `DataFrame` with one row for every unordered pair of sites and the
following columns:

- `site1`, `site2`: Identifiers of the two compared sites.
- `jaccard`: Jaccard dissimilarity, calculated from taxon presence/absence.
- `sorensen`: Sørensen dissimilarity, calculated from taxon presence/absence.
- `bray_curtis`: Bray–Curtis dissimilarity, calculated from the values in
  the community matrix.
- `baselga_total`: Total Sørensen dissimilarity used in the Baselga
  decomposition.
- `baselga_turnover`: Turnover component of the Baselga decomposition.
- `baselga_richness`: Richness-difference (nestedness-resultant) component
  of the Baselga decomposition.
- `carvalho_total`: Total Jaccard dissimilarity used in the Carvalho
  decomposition.
- `carvalho_turnover`: Replacement component of the Carvalho decomposition.
- `carvalho_richness`: Richness-difference component of the Carvalho
  decomposition.

# Calculation details

For each pair of sites, let:

- `a` be the number of taxa present in both sites;
- `b` be the number of taxa present only in the first site;
- `c` be the number of taxa present only in the second site.

A taxon is considered present when its matrix value is greater than zero.

The Jaccard and Sørensen dissimilarities are calculated as:

```text
Jaccard  = (b + c) / (a + b + c)
Sørensen = (b + c) / (2a + b + c)
```

Bray–Curtis is calculated from the numeric matrix values:

```text
Bray–Curtis = sum(abs.(x .- y)) / (sum(x) + sum(y))
```
where x and y are the value vectors for the two sites.

The Baselga decomposition partitions Sørensen dissimilarity into turnover
and richness-difference (nestedness-resultant) components:

```text
Baselga total     = (b + c) / (2a + b + c)
Baselga turnover  = min(b, c) / (a + min(b, c))
Baselga richness  = Baselga total - Baselga turnover
```

The Carvalho decomposition partitions Jaccard dissimilarity into turnover
(replacement) and richness-difference components:

```text
Carvalho total     = (b + c) / (a + b + c)
Carvalho turnover  = 2 * min(b, c) / (a + b + c)
Carvalho richness  = abs(b - c) / (a + b + c)
```

# Conventions and interpretation

Jaccard, Sørensen, Baselga and Carvalho calculations use presence/absence:
any value greater than zero is treated as presence. They do not account
for differences in abundance.

Bray–Curtis uses the matrix values directly. Its interpretation therefore
depends on the measurement semantics and scale represented by
c.values.

If both sites contain no recorded taxa, the corresponding dissimilarities
and components are returned as 0.0 to avoid division by zero. This is a
computational convention, not evidence that the sites are ecologically
identical.

A zero in the community matrix represents an unrecorded site–taxon value;
it does not necessarily indicate a verified ecological absence.
If fewer than two sites are present, the function returns an empty
DataFrame with the columns described above.

# References

- Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography, 19(1), 134–143. https://doi.org/10.1111/j.1466-8238.2009.00490.x

- Carvalho, J. C., Cardoso, P., Borges, P. A. V., Schmera, D., & Podani, J. (2013). Measuring fractions of beta diversity and their relationships to nestedness: A theoretical and empirical comparison of novel approaches. Oikos, 122(6), 825–834. https://doi.org/10.1111/j.1600-0706.2012.20980.x

# Example
results = beta_diversity(community)
"""
function beta_diversity(c::CommunityMatrix)
    # Initialise output dataframe
    d=DataFrame(
        site1=String[],             # Site ID 1
        site2=String[],             # Site ID 2
        jaccard=Float64[],          # Jaccard dissimilarity
        sorensen=Float64[],         # Sorensen dissimilarity
        bray_curtis=Float64[],      # Bray-curtis dissimilarity

        baselga_total=Float64[],    # Beta diversity according to Baselga framework
        baselga_turnover=Float64[],
        baselga_nestedness=Float64[],

        carvalho_total=Float64[],   # Beta diversity according to Carvalho framework
        carvalho_turnover=Float64[],
        carvalho_richness=Float64[]
    )

    nsites = length(c.sites)

    for i in 1:(nsites - 1), j in (i + 1):nsites
        a_values = Float64.(c.values[i, :])
        b_values = Float64.(c.values[j, :])

        # Presence/absence counts for the pair of sites
        shared = count(k -> a_values[k] > 0 && b_values[k] > 0,
                       eachindex(a_values))
        only_site1 = count(k -> a_values[k] > 0 && b_values[k] == 0,
                           eachindex(a_values))
        only_site2 = count(k -> a_values[k] == 0 && b_values[k] > 0,
                           eachindex(a_values))

        union_count = shared + only_site1 + only_site2
        total_abundance = sum(a_values) + sum(b_values)

        # Existing Jaccard and Sørensen dissimilarities
        jaccard = union_count == 0 ?
            0.0 : (only_site1 + only_site2) / union_count

        sorensen_denominator = 2 * shared + only_site1 + only_site2
        sorensen = sorensen_denominator == 0 ?
            0.0 : (only_site1 + only_site2) / sorensen_denominator

        # Existing abundance-based Bray–Curtis dissimilarity
        bray_curtis = total_abundance == 0 ?
            0.0 : sum(abs.(a_values .- b_values)) / total_abundance

        # Baselga: Sørensen dissimilarity partitioned into turnover
        # and nestedness-resultant (richness) components.
        baselga_total = sorensen
        baselga_turnover_denominator = shared + min(only_site1, only_site2)
        baselga_turnover = baselga_turnover_denominator == 0 ?
            0.0 : min(only_site1, only_site2) / baselga_turnover_denominator
        baselga_richness = baselga_total - baselga_turnover

        # Carvalho: Jaccard dissimilarity partitioned into replacement
        # and richness-difference components.
        carvalho_total = jaccard
        carvalho_turnover = union_count == 0 ?
            0.0 : 2 * min(only_site1, only_site2) / union_count
        carvalho_richness = union_count == 0 ?
            0.0 : abs(only_site1 - only_site2) / union_count

        push!(
            d,
            (
                c.sites[i],
                c.sites[j],
                jaccard,
                sorensen,
                bray_curtis,
                baselga_total,
                baselga_turnover,
                baselga_richness,
                carvalho_total,
                carvalho_turnover,
                carvalho_richness
            )
        )
    end

    return d
end

"""
    rank_abundance(c::CommunityMatrix)

Create a site-by-site rank-abundance table from a community matrix.

For each site, taxa are ranked by their recorded value in descending order.
The returned table contains one row for each taxon with a strictly positive
value. Its `proportion` is the taxon's value divided by the sum of values
recorded for that site.

This function is intended for quantitative abundance data. It throws an
`ArgumentError` unless `c.semantics == :abundance`.

# Arguments
- `c::CommunityMatrix`: A validated site-by-taxon matrix whose `semantics`
  must be `:abundance`. Rows represent sites and columns represent taxa.

# Returns
A `DataFrame` with the following columns:

- `site::String`: Site identifier.
- `rank::Int`: Rank of the taxon within its site, ordered from highest to
  lowest value. Ranks are assigned before zero-valued taxa are omitted.
- `taxon::String`: Taxon identifier.
- `value::Float64`: Recorded abundance for that site and taxon.
- `proportion::Float64`: The taxon's share of the total recorded abundance
  for that site (`value / sum(values)`).

Sites whose total abundance is zero are omitted. Taxa with zero values are
also omitted. For each included site, the proportions sum to approximately
one, subject to floating-point rounding.

# Errors
- Throws `ArgumentError` if `c.semantics` is not `:abundance`.

# Notes
- Tied values receive consecutive ranks; this function does not calculate
  shared or average ranks.
- The relative order of taxa with tied values is not explicitly defined.
- This function returns tabular data only; it does not create a plot.

# Example
```julia
results = rank_abundance(community)
```
"""
function rank_abundance(c::CommunityMatrix)
    # Ensure that community matrix is filled with abundance data
    c.semantics == :abundance ||
        throw(ArgumentError("rank_abundance requires a CommunityMatrix with semantics=:abundance"))
    # Initialise output dataframe
    d=DataFrame(
        site=String[],          # Site ID
        rank=Int[],             # Taxon rang in Site ID
        taxon=String[],         # Taxon name or identifier
        value=Float64[],        # Abundance value of the taxon in the site
        proportion=Float64[]    # Proportion of the taxon in the overall abundance of the site
    )

    # Computation are done by sites
    for (i,s) in enumerate(c.sites)
        x=c.values[i,:];                # Store abundance values
        n=sum(x);                       # Total abundance
        n==0&&continue                  # Ensure striclty positive abundance
        # Computation per taxon ordered following decreasing abundance
        for (r,j) in enumerate(sortperm(x;rev=true))
            # Add only taxa with striclty positive abundance value
            x[j]>0&&push!(
                d,
                (
                    s,
                    r,
                    c.taxa[j],
                    Float64(x[j]),
                    Float64(x[j]/n)
                )
            )
        end
    end
    d
end

"""
    accumulation_curve(c::CommunityMatrix; iterations=100, rng=Random.default_rng())

Estimate the species-accumulation curve by repeatedly adding sites in random
order. For each number of sites added, the function computes the cumulative
number of taxa observed and returns its mean and standard deviation across
iterations.

A taxon is considered observed at a site when its matrix value is greater
than zero. Thus, the function measures cumulative observed richness, regardless
of whether `c.semantics` is `:incidence`, `:abundance`, or `:records`.

# Arguments
- `c::CommunityMatrix`: Community matrix, with sites in rows and taxa in columns.
- `iterations::Integer=100`: Number of random site orderings to evaluate.
  Must be a positive integer.
- `rng=Random.default_rng()`: Random number generator used to shuffle site
  orderings. Supply a seeded generator, such as `MersenneTwister(1234)`, to
  make the random permutations reproducible.

# Returns
A `DataFrame` with one row for each number of sites added and these columns:

- `sites_added`: Number of sites included, from `1` to the number of sites.
- `mean_richness`: Mean cumulative richness across iterations.
- `sd_richness`: Standard deviation of cumulative richness across iterations.
  This is set to zero when `iterations == 1`.

# Notes
Each iteration starts with no taxa observed and visits every site exactly once,
in an order generated using `rng`. At each step, the function updates the set
of taxa observed so far and records its size. The returned statistics summarize
those values across iterations; they do not represent confidence intervals.

The input matrix is not modified. If the community contains no sites, the
returned `DataFrame` has empty columns.

# Throws
- `ArgumentError` if `iterations` is not a positive integer.

# Examples
```julia
using Random

curve = accumulation_curve(community; iterations=500)

# Use a seeded random number generator for reproducible permutations.
curve = accumulation_curve(
    community;
    iterations=500,
    rng=MersenneTwister(1234),
)
```
"""
function accumulation_curve(
    c::CommunityMatrix;
    iterations=100,
    rng=Random.default_rng() # random number generator used to shuffle the order of the sites in each repetition
    )
    # Check iterations number: need a striclty positive integer
    iterations isa Integer && iterations > 0 ||
        throw(ArgumentError("iterations must be a positive integer"))

    # Initialise output
    n = length(c.sites)             # Number of sites
    x = zeros(iterations, n)        # Cumulated richness

    # Repetition and sites shuffling
    for r in 1:iterations
        p = falses(length(c.taxa))
        for (k, i) in enumerate(randperm(rng, n))
            p .|= c.values[i, :] .> 0
            x[r, k] = count(p)
        end
    end

    DataFrame(
        sites_added=collect(1:n),
        mean_richness=n == 0 ? Float64[] : vec(mean(x, dims=1)),
        sd_richness=n == 0 ? Float64[] :
            iterations == 1 ? zeros(n) : vec(std(x, dims=1))
    )
end


"""
    richness_estimators(c::CommunityMatrix)

Estimate species richness using the estimator selected from the measurement
semantics of `c`:

- For `c.semantics == :abundance`, calculate Chao1 separately for each site.
- For `c.semantics == :incidence`, calculate Chao2 across all sites.

The function logs an informational message identifying the estimator used.
It does not modify `c`.

# Arguments
- `c::CommunityMatrix`: Community matrix with sites in rows and taxa in
  columns. Its `semantics` field determines which estimator is used.

# Returns
A `DataFrame`. Its columns depend on the selected estimator.

For `:abundance` semantics, the function returns one row per site with:

- `site`: Site identifier.
- `observed_richness`: Number of taxa with abundance greater than zero.
- `singletons`: Number of taxa represented by exactly one individual.
- `doubletons`: Number of taxa represented by exactly two individuals.
- `chao1`: Chao1 richness estimate for that site.
- `sample_coverage`: Coverage estimate, calculated as
  `1 - singletons / total_abundance`. It is `missing` when the site has no
  recorded abundance.

For `:incidence` semantics, the function returns one row for the whole
community with:

- `observed_richness`: Number of taxa recorded in at least one site.
- `uniques`: Number of taxa recorded in exactly one site.
- `duplicates`: Number of taxa recorded in exactly two sites.
- `sites_sampled`: Number of sites in the matrix.
- `chao2`: Chao2 richness estimate across the sites.

# Calculation
For abundance data, Chao1 is calculated per site as:

    observed_richness + singletons^2 / (2 * doubletons)

If there are no doubletons, the function uses this fallback correction:

    observed_richness + singletons * (singletons - 1) / 2

For incidence data, Chao2 is calculated across sites as:

    observed_richness + ((sites_sampled - 1) / sites_sampled)
        * uniques^2 / (2 * duplicates)

If there are no duplicates, the function uses this fallback correction:

    observed_richness + uniques * (uniques - 1) / (2 * sites_sampled)

# Requirements and interpretation
- Chao1 requires integer abundance counts. The function throws an
  `ArgumentError` if any value in an abundance matrix is not an integer.
- Chao2 uses presence or absence across sites. Positive values count as
  incidences, regardless of their magnitude.
- `:records` semantics are not supported and cause an `ArgumentError`.
- Chao2 is a community-wide estimate, not a separate estimate for each site.

# References
- Chao, A. (1984). Nonparametric estimation of the number of classes in a
  population. *Scandinavian Journal of Statistics*, 11(4), 265–270.
- Chao, A. (1987). Estimating the population size for capture-recapture data
  with unequal catchability. *Biometrics*, 43(4), 783–791.
- Colwell, R. K., & Coddington, J. A. (1994). Estimating terrestrial
  biodiversity through extrapolation. *Philosophical Transactions of the
  Royal Society of London. Series B: Biological Sciences*, 345(1311), 101–118.

# Throws
- `ArgumentError` if `c.semantics` is neither `:abundance` nor `:incidence`,
  or if abundance values are not integers.

# Examples
```julia
# Chao1: per-site estimates from integer abundance data
chao1_results = richness_estimators(abundance_community)

# Chao2: one community-wide estimate from incidence data
chao2_result = richness_estimators(incidence_community)
```
"""
function richness_estimators(c::CommunityMatrix)
    if c.semantics == :abundance

        # Accept floating-point abundances when they represent whole individuals.
        x_float = Float64.(c.values)

        all(x -> isfinite(x), x_float) ||
            throw(ArgumentError(
                "Chao1 requires finite abundances; fractional values cannot be converted."
            ))
        
        # Convert abundance values to integer
        x_integer = round.(Int, x_float)

        @info "Calculating Chao1 from abundance data."

        d = DataFrame(
            site=String[],
            observed_richness=Int[],
            singletons=Int[],
            doubletons=Int[],
            chao1=Float64[],
            sample_coverage=Union{Missing,Float64}[],
        )

        for (i, site) in enumerate(c.sites)
            # Convert validated integer-valued abundances before counting.
            x = x_integer[i, :]
            richness = count(>(0), x)
            singletons = count(==(1), x)
            doubletons = count(==(2), x)
            total = sum(x)

            chao1 = richness + (
                doubletons == 0 ?
                singletons * (singletons - 1) / 2 :
                singletons^2 / (2 * doubletons)
            )

            coverage = total == 0 ? missing : 1 - singletons / total

            push!(d, (site, richness, singletons, doubletons, chao1, coverage))
        end

        return d

    elseif c.semantics == :incidence
        @info "Calculating Chao2 from incidence data."

        n_sites = length(c.sites)
        site_frequencies = vec(sum(c.values .> 0, dims=1))
        observed_richness = count(>(0), site_frequencies)
        uniques = count(==(1), site_frequencies)
        duplicates = count(==(2), site_frequencies)

        chao2 = if n_sites == 0
            0.0
        elseif duplicates == 0
            observed_richness + uniques * (uniques - 1) / (2 * n_sites)
        else
            observed_richness + ((n_sites - 1) / n_sites) *
                (uniques^2 / (2 * duplicates))
        end

        return DataFrame(
            observed_richness=[observed_richness],
            uniques=[uniques],
            duplicates=[duplicates],
            sites_sampled=[n_sites],
            chao2=[Float64(chao2)],
        )

    else
        throw(ArgumentError(
            "Richness estimators require :abundance (Chao1) or :incidence (Chao2) semantics."
        ))
    end
end


"""
    sampling_coverage(c::CommunityMatrix)

Return the per-site sample coverage estimates.

The result contains only the `site` and `sample_coverage` columns from
[`richness_estimators`](@ref). The estimator and therefore the coverage
calculation depend on `c.semantics`: Chao1 is used for abundance data and
Chao2 for incidence data.

# Returns
A `DataFrame` with columns `site` and `sample_coverage`.
"""
function sampling_coverage(c::CommunityMatrix)
    estimates = richness_estimators(c)

    # Keep the output consistent if an estimator returns no coverage column.
    if !(:sample_coverage in propertynames(estimates))
        return DataFrame(
            site = copy(c.sites),
            sample_coverage = Union{Missing, Float64}[missing for _ in c.sites],
        )
    end

    return select(estimates, :site, :sample_coverage)
end


"""
Container for the tabular outputs of the supported community analyses.

Each field stores the result of one analysis as a `DataFrame`. `estimators`
is `nothing` when richness estimators are not applicable to the input
semantics.
"""
struct AnalysisResult
    alpha::DataFrame
    beta::DataFrame
    composition::DataFrame
    accumulation::DataFrame
    estimators::Union{Nothing, DataFrame}
    diagnostics::DataFrame
end


"""
    analyze(c::CommunityMatrix; iterations=100, rng=Random.default_rng(),
            diagnostics=DataFrame())

Run the supported analyses for a community matrix and collect their results.

Abundance matrices include a rank-abundance composition table and richness
estimators (Chao1). Incidence matrices do not run the rank-abundance analysis
and include richness estimators based on incidence (Chao2).

# Arguments
- `c`: Validated site-by-taxon community matrix.
- `iterations`: Number of random permutations used by the accumulation curve.
- `rng`: Random number generator passed to [`accumulation_curve`](@ref).
- `diagnostics`: Optional row-level diagnostics to retain in the result.

# Returns
An [`AnalysisResult`](@ref) containing the output tables.
"""
function analyze(
    c::CommunityMatrix;
    iterations = 100,
    rng = Random.default_rng(),
    diagnostics = DataFrame(),
)
    # Rank-abundance is defined only for abundance data.
    composition = if c.semantics == :abundance
        rank_abundance(c)
    else
        # Return an empty table with the same columns as rank_abundance.
        DataFrame(
            site = String[],
            rank = Int[],
            taxon = String[],
            value = Float64[],
            proportion = Float64[],
        )
    end

    # Richness estimators select Chao1 or Chao2 according to the semantics.
    estimators = c.semantics in (:abundance, :incidence) ?
        richness_estimators(c) :
        nothing

    return AnalysisResult(
        alpha_diversity(c),
        beta_diversity(c),
        composition,
        accumulation_curve(c; iterations, rng),
        estimators,
        DataFrame(diagnostics),
    )
end


"""
    export_tables(r::AnalysisResult, dir; overwrite=false)

Write the analysis result tables to CSV files in `dir`.

Files are written in a fixed order: `alpha.csv`, `beta.csv`,
`composition.csv`, `accumulation.csv`, optionally `estimators.csv`, and
`diagnostics.csv`. When `overwrite` is `false`, the function checks every
destination before writing and throws an `ArgumentError` if any already
exists.

# Returns
A vector of file paths in the same fixed order as the files are written.
"""
function export_tables(r::AnalysisResult, dir; overwrite = false)
    # An ordered vector of pairs guarantees deterministic write and return order.
    tables = Pair{String, DataFrame}[
        "alpha.csv" => r.alpha,
        "beta.csv" => r.beta,
        "composition.csv" => r.composition,
        "accumulation.csv" => r.accumulation,
    ]

    # Include estimators only when they were calculated.
    if r.estimators !== nothing
        push!(tables, "estimators.csv" => r.estimators)
    end

    push!(tables, "diagnostics.csv" => r.diagnostics)

    paths = [joinpath(dir, filename) for (filename, _) in tables]

    # Check all destinations before creating the directory or writing files.
    if !overwrite && any(ispath, paths)
        throw(ArgumentError("Refusing to overwrite existing tables"))
    end

    mkpath(dir)

    # Write files in the same stable order as `paths`.
    for (filename, table) in tables
        CSV.write(joinpath(dir, filename), table)
    end

    return paths
end