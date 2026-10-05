## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(ItemRest)

## ----simulate-----------------------------------------------------------------
set.seed(5102026)
n <- 600L
F1 <- rnorm(n)
F2 <- .3 * F1 + sqrt(1 - .3^2) * rnorm(n)
d <- data.frame(
  A1 = F1 + rnorm(n, sd=.5), A2 = F1 + rnorm(n, sd=.5),
  A3 = F1 + rnorm(n, sd=.5), A4 = F1 + rnorm(n, sd=.5),
  B1 = F2 + rnorm(n, sd=.5), B2 = F2 + rnorm(n, sd=.5),
  B3 = F2 + rnorm(n, sd=.5), B4 = -F2 + rnorm(n, sd=.5),
  Weak = rnorm(n), Cross = .6 * F1 + .6 * F2 + rnorm(n, sd=.5)
)
split <- itemrest_split(d, n_factors=2, seed=5102026,
  keys=c(B4=-1), retain_items=c("A1", "B1"),
  item_reasons=c(A1="Domain A anchor", B1="Domain B anchor"),
  max_solutions=500, verbose=FALSE)
result <- split$discovery
print(result)

## ----diagnostics--------------------------------------------------------------
result$removal_summary[, c("Solution_ID", "Baseline", "Removed_Items",
  "Analysis_Status", "Candidate", "Problem_Items", "Diagnostics")]
result$search
result$content_decisions

## ----reliability--------------------------------------------------------------
result$solution_details[["S00001"]]$assessment$factor_reliability
result$removal_summary[, c("Solution_ID", "Cronbachs_Alpha",
  "Standardized_Alpha", "Analysis_Correlation_Alpha", "Omega_Total")]

## ----holdout------------------------------------------------------------------
split$validation$validation_summary[, c("Solution_ID", "Candidate",
  "Analysis_Status", "Problem_Items", "N_Obs")]
split$validation$solution_details[["S00001"]]$factor_congruence

## ----cfa----------------------------------------------------------------------
if (requireNamespace("lavaan", quietly=TRUE)) {
  cfa <- itemrest_validate(result, d[split$validation_rows, ], method="cfa", seed=5102026)
  cfa$validation_summary
}

## ----bootstrap----------------------------------------------------------------
boot <- itemrest_bootstrap(result, d[split$discovery_rows, ], n_boot=20, seed=5102026)
boot$item_stability
boot$solution_stability
table(boot$replicates$Status)

