# Reproducible implementation example, not a claim of measurement superiority.
library(ItemRest)
set.seed(5102026)
n <- 600L
F1 <- rnorm(n)
F2 <- 0.3 * F1 + sqrt(1 - 0.3^2) * rnorm(n)
data <- data.frame(
  A1 = F1 + rnorm(n, sd = .5), A2 = F1 + rnorm(n, sd = .5),
  A3 = F1 + rnorm(n, sd = .5), A4 = F1 + rnorm(n, sd = .5),
  B1 = F2 + rnorm(n, sd = .5), B2 = F2 + rnorm(n, sd = .5),
  B3 = F2 + rnorm(n, sd = .5), B4 = -F2 + rnorm(n, sd = .5),
  Weak = rnorm(n), Cross = .6 * F1 + .6 * F2 + rnorm(n, sd = .5)
)
split <- itemrest_split(
  data, train_fraction = .7, seed = 5102026, n_factors = 2,
  keys = c(B4 = -1), retain_items = c("A1", "B1"),
  item_reasons = c(A1 = "Domain A anchor", B1 = "Domain B anchor"),
  max_solutions = 500, rank_by = "none", verbose = FALSE
)
print(split$discovery)
print(split$discovery, report = "all")
split$validation$validation_summary
split$discovery$solution_details[["S00001"]]$assessment$factor_reliability
if (requireNamespace("lavaan", quietly = TRUE)) {
  cfa <- itemrest_validate(split$discovery, data[split$validation_rows, ], method = "cfa", seed = 5102026)
  print(cfa$validation_summary)
}
bootstrap <- itemrest_bootstrap(split$discovery, data[split$discovery_rows, ], n_boot = 20, seed = 5102026)
bootstrap$item_stability
bootstrap$replicates
split$discovery$provenance
