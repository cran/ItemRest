# ItemRest

ItemRest 1.0.0 provides a threshold-driven search for candidate item sets in
exploratory factor analysis. It documents removal paths and numerical problems
and supports independent evaluation. Candidate eligibility does not establish
an optimal, valid, or stable measurement instrument.

## Installation

Install this updated local source tree from a terminal:

```sh
R CMD INSTALL path/to/ItemRest
```

Repository: [ItemRest](https://github.com/ahmetcaliskan1987/ItemRest).
Published releases may precede this local update. License: MIT.

## Candidate search

```r
library(ItemRest)
result <- itemrest(
  data = my_data,
  n_factors = 3,
  missing = "listwise",
  keys = c(reverse_item = -1),
  retain_items = "content_anchor",
  item_reasons = c(content_anchor = "Required domain coverage"),
  rank_by = "none"
)
result$candidate_solutions
result$removal_summary
print(result, report = "all")
```

Replace the example item names with names in your data, or omit these arguments.
Each removal is followed by reassessment of the remaining items. Only flagged
items are removable; this is not a search over every possible subset of all
items. Duplicate retained sets are evaluated once. Factor count is fixed at
baseline, including when parallel analysis determines it. The baseline always
appears first in the full summary, even if it is not a candidate.

The correlation matrix and pairwise observation counts are computed once across
all baseline items, after missing-data handling and scoring keys. Parallel
analysis and each removal step use the appropriate submatrix. Positive
definiteness and requested smoothing are checked separately for each retained
set, using the original unsmoothed correlations. The original matrices are
available as `result$correlation_matrix` and `result$pairwise_n`. Validation EFA
computes its own matrix, and each bootstrap replicate computes a new matrix.

Extraction and rotation names are checked before analysis. Supported mixed-case
names such as `Oblimin` are canonicalized. Use `geominQ` or `geominT` explicitly;
ambiguous or unsupported names raise an error. The default oblimin rotation's
`GPArotation` dependency is declared. Estimator messages are saved in `Messages`.
Informational messages do not change screening eligibility. Explicit reports of
nonconvergence, unavailable rotations, variance problems, or matrix repair
require review. Source-matrix messages and warnings are recorded once as
`correlation_messages` and `correlation_warnings`, not repeated for each subset.
The ordinal category limit is configurable through `ordinal_categories = 7`;
missing values are excluded from the count.

`parallel_method = "fa"` remains the default, using reduced-matrix eigenvalues
from a one-factor minres fit. Choose `parallel_method = "pc"` for full-matrix
component eigenvalues. For polychoric/mixed correlations, reference samples
independently permute each item's observed values, preserving category margins
and missing positions, and use the same correlation estimator. This permutation
reference differs from simulating normal Pearson data and can change the suggested
factor count. Inspect `result$parallel_analysis`; supplying `n_factors` bypasses PA.
Computing reference matrices costs time; the real-data matrix is still computed
only once for the removal search.

Loadings, Phi, communalities, reliability, and diagnostics are retained by default.
Use `store_fits = TRUE` to also retain full psych EFA fit objects. Protecting an
item with `retain_items` does not exempt it from screening; a flagged protected
item can prevent any candidate solutions.

Validation EFA matches factor order and signs by maximizing total absolute
congruence and reports `Min_Congruence`. CFA also reports scaled fit indices.
Identical-data checks compare numeric values and item names, independent of row
names and integer/double storage; they cannot prove independent observations.
Validation also detects reordered copies of discovery responses, while bootstrap
requires the original row order. Validation EFA correlates the union of requested
retained items. Constant items fail the sets containing them; other requested
sets can still be evaluated from one matrix of the varying items. Missing-data
preparation continues to fix rows across all discovery items.
Bootstrap adds `*_All_Attempts` frequencies with all attempted replicates in
the denominator, treating failed/incomplete runs as having no candidates. These
are conservative sensitivity summaries alongside the conditional frequencies.
`progress = TRUE` displays a bootstrap progress bar; it is enabled interactively.

The default seed is the local execution date in DDMMYYYY form: 5 October 2026
produces label `05102026` and numeric seed `5102026`. Supply an explicit seed
for reproduction on another date. Seeds, thresholds, scoring keys, R/package
versions, retained observations, and search completeness are recorded. Supplied
and automatic seeds restore the caller's RNG state. Parallel analysis uses
100 replications by default, controlled by `parallel_iterations`.

No single winner is selected. The default order is discovery order. Optional
`rank_by = "n_removed"` or `"explained_variance"` changes presentation while
keeping the full-summary baseline first. Explained variance is mean model
communality for the retained subset; differences across subsets do not establish
measurement superiority. `max_solutions` defaults to 10000, and limited searches
are explicitly marked incomplete.

## Diagnostics and reliability

Each attempted solution records status, errors, warnings, original correlation
positive definiteness, smoothing, available convergence information, communalities,
uniquenesses, factor correlations, factor support, and sample sizes. Failed or
unidentified branches do not abort the other branches. Convergence is labelled
`Not_reported` when the estimator exposes no convergence flag and no convergence
warning, rather than claiming convergence. A smoothed matrix, large pattern
coefficients, high factor correlations, or estimation warnings require review.

By default candidates need three qualifying primary items per factor, no flagged
loadings, numerical admissibility, and no review flags. These are configurable
screening safeguards. A loading above one requires inspection, especially with
oblique rotation; it is not equated automatically with a Heywood case.

`solution_details[["S00001"]]$assessment$factor_reliability` provides factor-level
results. Factor reliability uses items with qualifying, nonflagged primary
assignments. Global raw alpha, standardized Pearson alpha, selected-correlation
alpha, and model omega total are separately labelled. Raw alpha remains available
as `Cronbachs_Alpha` for table compatibility; it is not ordinal alpha. Model omega
total uses the fitted common covariance and uniquenesses of a standardized
unit-weighted sum. It includes all common factors and does not establish
unidimensionality. No automatic reverse scoring occurs: supply `keys` explicitly.

Listwise deletion is performed once across all baseline items. Pairwise handling
records the minimum/maximum pair counts and uses the minimum as conservative
EFA sample size. Neither option guarantees that missing-data bias is resolved.

## Independent evaluation and bootstrap stability

```r
validation <- itemrest_validate(result, independent_data)
validation$validation_summary

# Optional lavaan CFA: fixed discovery primary-factor assignments.
cfa <- itemrest_validate(result, independent_data, method = "cfa")

# Alternatively, split observations before searching.
split <- itemrest_split(my_data, n_factors = 3, seed = 5102026)

# Repeat the entire search, not just a previously selected solution.
stability <- itemrest_bootstrap(result, my_data, n_boot = 100, seed = 5102026)
stability$item_stability
stability$replicates
```

Validation never removes additional items. EFA includes factor-congruence matrices
for inspection. CFA fixes discovery primary assignments and zero cross-loadings,
uses MLR for continuous outcomes or WLSMV for specified/detected ordinal outcomes,
and reports fit indices without a universal pass/fail rule. Continuous CFA uses
FIML when the discovery missing policy is pairwise; this difference is recorded.
The user must supply independent observations. Refitting discovery data is not
independent validation.

Bootstrap holds the discovery factor count fixed and reruns the full search.
Item frequencies distinguish retention in any candidate, all candidates, and
average retention among candidates. Exact-set frequencies avoid rank-based
winner claims. Failed and incomplete replications are reported and excluded from
frequency denominators. Stability frequencies are not probabilities of validity.

## Reproducible example and migration

See the package vignette and `inst/examples/reproducible_workflow.R` for a complete
simulated discovery/holdout/bootstrap workflow. It is an implementation example,
not empirical evidence that this procedure improves measurement validity.

Version 1.0.0 replaces `optimal_strategy` with the complete
`candidate_solutions` table. `print(result, report = "candidates")` is the default.
The old `report = "optimal"` value warns and redirects to the candidate report.
Numeric item names are preserved, including periods and spaces. There is no
Shiny dependency or application.
