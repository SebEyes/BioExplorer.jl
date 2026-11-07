# Methods and interpretation

Richness counts strictly positive taxa. Shannon uses natural logarithms. Simpson is concentration `sum(p^2)`, not its complement. Hill q=1 is `exp(Shannon)`, q=2 is inverse concentration; both are zero for empty sites. Pielou is missing for richness below two.

Jaccard and Sørensen operate on positive-value support. Bray–Curtis uses the actual values and therefore inherits their semantics. Empty–empty dissimilarity is defined as zero. These conventions do not turn unrecorded taxa into verified absences.

Accumulation independently permutes sites for each replicate without changing the input. Supply an RNG for reproducibility. Reported SD describes random-order variability, not a confidence interval. One iteration returns zero SD.

Chao1 requires abundance semantics and integer counts. The estimator is S + f1²/(2f2) when f2 > 0, otherwise S + f1(f1−1)/2. Good's singleton coverage is 1 − f1/N, missing when N=0; it estimates observed probability mass, not the fraction of all species observed. No coverage-standardized rarefaction is implemented. See Chao (1984), Scandinavian Journal of Statistics 11:265–270, and Good (1953), Biometrika 40:237–264.

Darwin Core ingestion does not resolve taxonomy, deduplicate occurrence IDs, infer sampling effort, or parse full DwC-A archives. Explicit absence rows are excluded and cannot establish an absence-only site's sampling design. Missing numeric counts are never imputed. Users must justify grouping, comparability, and the ecological interpretation of zeros.
