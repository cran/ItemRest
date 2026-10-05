# ItemRest 1.0.0

* Record informational estimator messages without disqualifying candidates.
  Silently preload GPArotation for rotations and restrict message-driven review
  to explicit convergence, rotation, variance, and matrix-repair problems.
* Scope validation correlations to requested item sets and isolate constant-item
  failures. Suppress diagnostic variable tables printed on correlation errors.
* Record source correlation warnings/messages once per analysis, rather than
  attributing them to every subset. Detect identical validation response values
  after row reordering; retain order-sensitive bootstrap input checks.
* Label unassessed whole-set alpha as Not_assessed and remove unused combination
  generation, sorting helpers, and the obsolete internal combs argument.
* Validate extraction and rotation names before analysis, canonicalize supported
  capitalization, and declare the default rotation's GPArotation dependency.
  Record estimator messages and flag reported nonconvergence.
* Label Howard primary loadings below .40 as low loadings, and match holdout
  factors and signs using linear assignment rather than factorial enumeration.
  Add Min_Congruence and scaled CFA fit indices.
* Fingerprint numeric values and item names independently of row names and
  integer/double storage. Explain that protected flagged items can prevent
  candidate solutions and that unknown convergence remains Not_reported.
* Keep fa parallel analysis as default and add parallel_method = "pc". Ordinal
  reference samples preserve observed category margins and missing positions
  and use the same mixed/polychoric correlation estimator as the observed data.
  Add configurable ordinal_categories, excluding NA from category counts.
* Add bootstrap frequencies across all attempts as a conservative sensitivity
  summary, predraw replicate rows and seeds, and provide a progress bar.
* Label whole-set alpha as descriptive in multifactor models; factor-specific
  reliability remains available. Store compact EFA results by default, with
  store_fits = TRUE for full psych fit objects.
* Compute the original correlation matrix and pairwise observation counts once
  per dataset and reuse submatrices in parallel analysis and removal searches.
  Check positive definiteness and apply requested smoothing per retained set.
  Validation EFA uses its own matrix; each bootstrap replicate uses a new matrix.
  Expose correlation_matrix and pairwise_n in discovery results.
* Replace the single optimal_strategy field with all screened candidate_solutions.
  Candidate reporting is the default; report = "optimal" is a deprecated alias.
* Reassess remaining items after removal, cache identical item sets, retain the
  no-removal baseline, and label bounded searches as incomplete.
* Add numerical admissibility, correlation positive-definiteness and explicit
  smoothing, convergence status, factor correlations, factor support, warnings,
  and branch-specific failure records. Unidentified subsets are skipped.
* Add explicit missing-data handling, effective pair counts, parameter checks,
  protected items, content notes, and scoring keys.
* Default ordering is discovery order; optional removal-count/variance ordering
  never selects a winner. Oblique explained variance uses mean model communality.
* Add factor-level alpha and model omega total, with separate raw, standardized,
  and selected-correlation alpha labels. Preserve signed loading ranges and add
  absolute ranges. Capture Howard failures below the secondary .30 threshold.
* Add holdout EFA, optional lavaan CFA, disjoint sample splitting, and bootstrap
  stability of the entire search with explicit failure/incompleteness denominators.
* Add DDMMYYYY date-based default seeds (e.g. label 05102026 / seed 5102026),
  caller RNG restoration, configurable parallel-analysis replications, and
  recorded settings and software provenance.
* Provide a reproducible simulated vignette and complete executable workflow.
* Restore "Let algorithms be your compass, not your captain." as the final line
  of every printed analysis summary.

# ItemRest 0.2.5.9000 (development version)

* Reassess low-loading and cross-loading items after each removal and test
  further combinations of the remaining flagged items. Each distinct item set
  is evaluated once and the factor count stays fixed throughout the search.
* Add iteration, step-specific and cumulative removal counts, remaining items,
  and remaining problem items to the removal summary while retaining existing
  columns.
* Require both low-loading and cross-loading problems to be resolved for optimal
  strategy selection and reporting.

# ItemRest 0.2.5

* Fixed: The factor loading range in the summary table now correctly displays raw values (including negative signs) instead of absolute values.
