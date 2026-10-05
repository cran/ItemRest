test_that("informational messages do not change candidate eligibility", {
  data <- fixture_data()
  count <- 0L
  local_mocked_bindings(efa_custom = function(data, ...) {
    count <<- count + 1L
    out <- make_test_efa(names(data))
    out$fit_messages <- if (count == 1L)
      c("Loading required namespace: GPArotation", "Informational estimator message") else character()
    out
  })
  first <- itemrest(data, n_factors = 2, rotate = "promax", reliability = "none", verbose = FALSE)
  second <- itemrest(data, n_factors = 2, rotate = "promax", reliability = "none", verbose = FALSE)
  expect_true(first$removal_summary$Candidate)
  expect_identical(first$candidate_solutions$Solution_ID, second$candidate_solutions$Solution_ID)
  expect_false(first$removal_summary$Review_Required)
  expect_equal(first$removal_summary$Diagnostics, "")
  expect_match(first$removal_summary$Messages, "Informational estimator message")
  expect_length(estimation_message_flags("Loading required namespace: GPArotation"), 0L)
  expect_equal(estimation_message_flags("Specified rotation not found, rotate='none' used"), "rotation_unavailable")
  expect_equal(estimation_message_flags("A Heywood case was detected"), "reported_variance_problem")
})

test_that("actual promax calls give the same candidates and numerical results", {
  data <- with_itemrest_seed(29, {
    factors <- matrix(rnorm(800), ncol = 2)
    as.data.frame(sapply(rep(1:2, each = 4), function(i) factors[, i] + rnorm(400, sd = .4)))
  })
  first <- itemrest(data, n_factors = 2, rotate = "promax", reliability = "none", verbose = FALSE, seed = 17)
  second <- itemrest(data, n_factors = 2, rotate = "promax", reliability = "none", verbose = FALSE, seed = 17)
  expect_equal(nrow(first$candidate_solutions), 1L)
  expect_equal(first$removal_summary, second$removal_summary)
  expect_equal(first$initial_efa$loadings, second$initial_efa$loadings)
})

test_that("validation only correlates selected items and isolates constant-item failures", {
  data <- real_factor_data(n = 400, p = 7)
  data <- as.data.frame(lapply(data, function(x) as.integer(cut(x, c(-Inf, -.7, 0, .7, Inf)))))
  local_mocked_bindings(identify_problem_items = function(efa_res, ...) {
    problems <- empty_problem_items()
    problems$low_loading <- intersect("V1", rownames(efa_res$loadings))
    problems
  })
  discovery <- itemrest(data[1:200, ], cor_method = "polychoric", n_factors = 1,
    reliability = "none", verbose = FALSE)
  candidate <- discovery$candidate_solutions$Solution_ID
  expect_length(candidate, 1L)
  validation <- data[201:400, ]
  validation$V1 <- 1L
  inputs <- list()
  original_cor <- cor_matrix_custom
  local_mocked_bindings(cor_matrix_custom = function(data, ...) {
    inputs[[length(inputs) + 1L]] <<- data
    original_cor(data, ...)
  })
  selected <- itemrest_validate(discovery, validation, solution_ids = candidate)
  expect_equal(selected$validation_summary$Analysis_Status, "ok")
  expect_true(selected$validation_summary$Candidate)
  expect_length(inputs, 1L)
  expect_false("V1" %in% names(inputs[[1]]))
  expect_equal(selected$correlation_items, names(validation)[-1L])
  inputs <- list()
  all <- itemrest_validate(discovery, validation)
  expect_equal(all$validation_summary$Analysis_Status, c("failed", "ok"))
  expect_match(all$validation_summary$Error_Message[1], "no variation")
  expect_equal(all$constant_items, "V1")
  expect_equal(all$validation_summary$Global_Alpha_Interpretation[1], "Not_assessed")
  expect_length(inputs, 1L)
})

test_that("source correlation warnings are recorded once without attributing them to subsets", {
  data <- real_factor_data(n = 150, p = 7)
  original_cor <- cor_matrix_custom
  calls <- 0L
  local_mocked_bindings(cor_matrix_custom = function(...) {
    calls <<- calls + 1L
    warning("Warning about a pair in the original correlation matrix")
    message("Original correlation information")
    original_cor(...)
  }, identify_problem_items = function(efa_res, ...) {
    problems <- empty_problem_items()
    problems$low_loading <- intersect("V1", rownames(efa_res$loadings))
    problems
  })
  result <- itemrest(data, n_factors = 1, reliability = "none", verbose = FALSE)
  expect_equal(calls, 1L)
  expect_equal(result$correlation_warnings, "Warning about a pair in the original correlation matrix")
  expect_equal(result$correlation_messages, "Original correlation information")
  expect_equal(result$removal_summary$Warnings, rep("", 2L))
  expect_equal(result$removal_summary$Messages, rep("", 2L))
  expect_true(result$removal_summary$Candidate[2L])
})

test_that("validation detects reordered response multisets while bootstrap requires row order", {
  data <- real_factor_data(n = 150)
  discovery <- itemrest(data, n_factors = 1, reliability = "none", verbose = FALSE)
  reordered <- data[rev(seq_len(nrow(data))), ]
  expect_warning(check <- itemrest_validate(discovery, reordered), "not independent")
  expect_equal(check$independence_note, "same_as_discovery")
  expect_error(itemrest_bootstrap(discovery, reordered, n_boot = 2), "must match")
  duplicate <- data[c(1:149, 1), ]
  expect_false(identical(data_fingerprint(duplicate, TRUE), data_fingerprint(data, TRUE)))
  values <- data.frame(A = c(2, NA, 1, 1), B = c(NA, 3, 4, 4))
  expect_identical(data_fingerprint(values, TRUE), data_fingerprint(values[c(4, 2, 1, 3), ], TRUE))
})

test_that("lavaan variable tables do not leak from failed correlation calls", {
  data <- data.frame(V1 = rep(1L, 50), V2 = rep(1:5, 10), V3 = rep(1:2, 25))
  output <- capture.output(result <- capture_analysis(cor_matrix_custom(data, "polychoric", "pairwise")))
  expect_length(output, 0L)
  expect_false(is.null(result$error))
})
