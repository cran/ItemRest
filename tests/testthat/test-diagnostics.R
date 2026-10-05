test_that("numerical and support problems prevent candidate eligibility", {
  data <- fixture_data()
  opts <- analysis_options("listwise", 3L, 0.85, "fail", "none")
  evaluate <- function(out) {
    local_mocked_bindings(efa_custom = function(...) out)
    evaluate_item_set(data, 2, "pearson", "uls", "oblimin", 0.3, 0.1, opts)
  }
  heywood <- make_test_efa(names(data))
  heywood$loadings[1, 1] <- 1.10
  result <- evaluate(heywood)
  expect_equal(result$status, "inadmissible")
  expect_false(result$assessment$candidate)
  expect_true("inadmissible_variance" %in% result$assessment$notes)
  high_phi <- evaluate(make_test_efa(names(data), phi = matrix(c(1, .9, .9, 1), 2)))
  expect_true(high_phi$assessment$review)
  expect_false(high_phi$assessment$candidate)
  failed <- make_test_efa(names(data))
  failed$fit_warnings <- "failed to converge"
  expect_equal(evaluate(failed)$assessment$convergence, "No")
  support <- make_test_efa(names(data))
  support$loadings[, 1] <- .7
  support$loadings[, 2] <- .05
  expect_false(evaluate(support)$assessment$candidate)
})

test_that("nonpositive correlations fail or are explicitly smoothed", {
  data <- fixture_data(p = 4)
  bad <- matrix(.9, 4, 4); diag(bad) <- 1; bad[1, 2] <- bad[2, 1] <- -.9
  local_mocked_bindings(cor_matrix_custom = function(...) bad)
  failed <- itemrest(data, n_factors = 1, reliability = "none", verbose = FALSE)
  expect_equal(failed$removal_summary$Analysis_Status, "failed")
  expect_match(failed$removal_summary$Error_Message, "positive definite")
  expect_false(failed$removal_summary$Positive_Definite)
  expect_lt(failed$removal_summary$Min_Eigenvalue, 0)
  smooth <- itemrest(data, n_factors = 1, reliability = "none", verbose = FALSE, pd_action = "smooth")
  expect_true(smooth$removal_summary$Correlation_Adjusted[1L])
  expect_true(smooth$removal_summary$Review_Required[1L])
  expect_false(smooth$removal_summary$Candidate[1L])
})

test_that("automatic factor selection is reproducible and preserves RNG", {
  data <- real_factor_data(n = 180, p = 6)
  set.seed(246)
  before <- .Random.seed
  first <- itemrest(data, verbose = FALSE, seed = 81, parallel_iterations = 20)
  second <- itemrest(data, verbose = FALSE, seed = 81, parallel_iterations = 20)
  expect_equal(first$settings$n_factors, second$settings$n_factors)
  expect_equal(first$removal_summary, second$removal_summary)
  expect_identical(.Random.seed, before)
  expect_equal(first$settings$factor_count_policy, "fixed_at_baseline")
})

test_that("seeded code restores an initially absent RNG state", {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (existed) rm(list = ".Random.seed", envir = .GlobalEnv)
  with_itemrest_seed(12, runif(3))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("missing policy records effective sample sizes and fixes listwise rows", {
  data <- real_factor_data(n = 120)
  data[1:5, 1] <- NA; data[6:10, 2] <- NA
  expect_error(itemrest(data, n_factors = 1, missing = "fail"), "Missing")
  listwise <- itemrest(data, n_factors = 1, verbose = FALSE)
  pairwise <- itemrest(data, n_factors = 1, missing = "pairwise", verbose = FALSE)
  expect_equal(listwise$data_summary$n_used, 110L)
  expect_equal(listwise$removal_summary$N_Obs, 110)
  expect_equal(pairwise$data_summary$n_used, 120L)
  expect_equal(pairwise$removal_summary$N_Obs, 110)
  expect_equal(pairwise$removal_summary$N_Obs_Max, 120)
})

test_that("factor reliability and omega distinguish scoring direction", {
  data <- real_factor_data()
  data$V1 <- -data$V1
  unkeyed <- itemrest(data, n_factors = 1, verbose = FALSE)
  keyed <- itemrest(data, n_factors = 1, keys = c(V1 = -1), verbose = FALSE)
  expect_gt(keyed$removal_summary$Cronbachs_Alpha, unkeyed$removal_summary$Cronbachs_Alpha)
  expect_equal(keyed$settings$keys["V1"], c(V1 = -1))
  factor <- keyed$solution_details[[1]]$assessment$factor_reliability
  expect_equal(factor$N_Items, 5L)
  expect_true(is.finite(factor$Omega_Total))
  expect_gt(factor$Omega_Total, 0)
  expect_lt(factor$Omega_Total, 1)
  expect_match(unkeyed$removal_summary$Factor_Loading_Range, "-")
  expect_false(grepl("-", unkeyed$removal_summary$Absolute_Loading_Range))
  L <- matrix(c(.6, .7, .8), ncol = 1)
  expected <- sum(L)^2 / (sum(L)^2 + sum(1 - L^2))
  expect_equal(model_omega(L, matrix(1), 1 - L^2), expected)
})

test_that("seeds use DDMMYYYY by default and preserve caller RNG", {
  expected <- format(Sys.Date(), "%d%m%Y")
  info <- resolve_itemrest_seed(NULL)
  expect_equal(info$label, expected)
  expect_equal(info$value, as.integer(expected))
  set.seed(93)
  before <- .Random.seed
  data <- real_factor_data()
  first <- itemrest(data, n_factors = 1, verbose = FALSE)
  second <- itemrest(data, n_factors = 1, verbose = FALSE)
  expect_identical(.Random.seed, before)
  expect_equal(first$removal_summary, second$removal_summary)
  expect_equal(first$settings$seed_label, expected)
  expect_equal(itemrest(data, n_factors = 1, seed = 12, verbose = FALSE)$settings$seed, 12L)
  expect_error(itemrest(data, seed = -1), "seed")
  expect_error(itemrest(data, n_factors = 0), "n_factors")
  expect_error(itemrest(data, min_loading = 2), "min_loading")
  expect_error(itemrest(data, keys = c(V1 = 0)), "keys")
})
