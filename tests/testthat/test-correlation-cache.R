test_that("removal searches reuse one matrix and match independent refits", {
  continuous <- real_factor_data(n = 220, p = 7)
  original_cor <- cor_matrix_custom
  local_mocked_bindings(identify_problem_items = function(efa_res, ...) {
    problems <- empty_problem_items()
    problems$low_loading <- intersect(c("V1", "V2"), rownames(efa_res$loadings))
    problems
  })
  for (method in c("pearson", "spearman", "kendall", "polychoric")) {
    data <- if (method == "polychoric") as.data.frame(lapply(continuous, function(x)
      as.integer(cut(x, c(-Inf, -.7, 0, .7, Inf))))) else continuous
    data[1:8, 1] <- NA
    data[9:14, 2] <- NA
    for (missing in c("listwise", "pairwise")) {
      calls <- 0L
      local_mocked_bindings(cor_matrix_custom = function(...) {
        calls <<- calls + 1L
        original_cor(...)
      })
      result <- itemrest(data, cor_method = method, n_factors = 1, missing = missing,
        keys = c(V3 = -1), reliability = "none", seed = 13, verbose = FALSE)
      expect_equal(calls, 1L)
      expect_equal(nrow(result$removal_summary), 4L)
      prepared <- apply_item_keys(prepare_itemrest_data(data, missing)$data, c(V3 = -1))$data
      expect_equal(result$correlation_matrix, as.matrix(original_cor(prepared, method, missing)))
      expect_equal(result$pairwise_n, pairwise_counts(prepared))
      for (detail in result$solution_details) {
        items <- detail$remaining_items
        independent <- efa_custom(prepared[, items, drop = FALSE], n_factors = 1,
          cor_method = method, missing = missing, reliability = "none")
        expect_equal(detail$analysis$correlation$matrix, independent$correlation$matrix,
          tolerance = 1e-7)
        expect_equal(detail$analysis$loadings, independent$loadings, tolerance = 1e-7)
        expect_equal(detail$analysis$explained_var, independent$explained_var, tolerance = 1e-7)
        expect_equal(detail$analysis$pairwise_n, independent$pairwise_n)
      }
    }
  }
})

test_that("parallel analysis and EFA share the baseline calculation", {
  calls <- 0L
  original_cor <- cor_matrix_custom
  local_mocked_bindings(cor_matrix_custom = function(...) {
    calls <<- calls + 1L
    original_cor(...)
  })
  result <- itemrest(real_factor_data(n = 200, p = 6), parallel_iterations = 10,
    seed = 32, reliability = "none", verbose = FALSE)
  expect_equal(calls, 1L)
  expect_equal(result$settings$correlation_policy, "baseline_submatrix")
  expect_equal(result$initial_efa$correlation$matrix, result$correlation_matrix)
})

test_that("smoothing is local to subsets and does not change the cached matrix", {
  data <- fixture_data(p = 4)
  raw <- matrix(.9, 4, 4, dimnames = list(names(data), names(data)))
  diag(raw) <- 1
  raw[1, 2] <- raw[2, 1] <- -.9
  calls <- 0L
  local_mocked_bindings(cor_matrix_custom = function(...) {
    calls <<- calls + 1L
    raw
  })
  cache <- new_correlation_cache(data, "pearson", "listwise")
  expect_warning(baseline <- checked_correlation(data, "pearson", "listwise", "smooth", cache),
    "smoothing was done")
  expect_true(baseline$adjusted)
  items <- names(data)[-1L]
  subset <- checked_correlation(data[, items], "pearson", "listwise", "smooth", cache)
  expect_false(subset$adjusted)
  expect_true(subset$positive_definite)
  expect_equal(subset$matrix, raw[items, items])
  expect_equal(cached_correlation_matrix(cache), raw)
  expect_equal(calls, 1L)
})

test_that("cached failures and warnings remain visible without recomputation", {
  data <- fixture_data()
  calls <- 0L
  local_mocked_bindings(cor_matrix_custom = function(...) {
    calls <<- calls + 1L
    warning("warning before correlation failure")
    stop("correlation estimation failed")
  })
  cache <- new_correlation_cache(data, "pearson", "listwise")
  opts <- analysis_options("listwise", 3L, .85, "fail", "none")
  evaluated <- lapply(list(data, data[, -1], data[, -2]), function(subset)
    evaluate_item_set(subset, 2, "pearson", "uls", "oblimin", .3, .1, opts, cache))
  expect_equal(calls, 1L)
  expect_equal(vapply(evaluated, `[[`, character(1), "status"), rep("failed", 3L))
  expect_true(all(vapply(evaluated, function(x) grepl("correlation estimation failed", x$error), logical(1))))
  expect_true(all(vapply(evaluated, function(x) length(x$warnings) == 0L, logical(1))))
  expect_equal(cache$result$warnings, "warning before correlation failure")
  local_mocked_bindings(cor_matrix_custom = function(...) {
    calls <<- calls + 1L
    warning("source correlation warning")
    stats::cor(data)
  })
  calls <- 0L
  cache <- new_correlation_cache(data, "pearson", "listwise")
  baseline <- capture_analysis(checked_correlation(data, "pearson", "listwise", "fail", cache))
  subset <- capture_analysis(checked_correlation(data[, -1], "pearson", "listwise", "fail", cache))
  expect_length(baseline$warnings, 0L)
  expect_length(subset$warnings, 0L)
  expect_equal(cache$result$warnings, "source correlation warning")
  expect_equal(calls, 1L)
})

test_that("validation and bootstrap calculate new matrices for their own observations", {
  data <- real_factor_data(n = 240, p = 7)
  original_cor <- cor_matrix_custom
  calls <- 0L
  local_mocked_bindings(cor_matrix_custom = function(...) {
    calls <<- calls + 1L
    original_cor(...)
  }, identify_problem_items = function(efa_res, ...) {
    problems <- empty_problem_items()
    problems$low_loading <- intersect(c("V1", "V2"), rownames(efa_res$loadings))
    problems
  })
  training <- itemrest(data[1:120, ], n_factors = 1, reliability = "none", verbose = FALSE)
  expect_equal(calls, 1L)
  ids <- names(training$solution_details)
  validation <- itemrest_validate(training, data[121:240, ], solution_ids = ids)
  expect_equal(calls, 2L)
  expect_equal(nrow(validation$validation_summary), 4L)
  expect_equal(validation$solution_details[[1]]$out$correlation$matrix, stats::cor(data[121:240, ]))
  bootstrap <- itemrest_bootstrap(training, data[1:120, ], n_boot = 3, seed = 17)
  expect_equal(calls, 5L)
  expect_equal(bootstrap$settings$complete_replicates, 3L)
  expect_equal(bootstrap$replicates$N_Evaluated, rep(4L, 3L))
})
