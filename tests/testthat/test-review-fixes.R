test_that("unsupported methods fail early and mixed-case names are canonicalized", {
  data <- fixture_data()
  expect_error(itemrest(data, extract = "xyz", n_factors = 2), "Unsupported extract")
  expect_error(itemrest(data, rotate = "geomin", n_factors = 2), "Unsupported rotate")
  methods <- validate_efa_methods("ULS", "Oblimin")
  expect_equal(methods, list(extract = "uls", rotate = "oblimin"))
  expect_equal(validate_efa_methods("ml", "Promax")$rotate, "Promax")
  expect_error(itemrest(data, n_factors = 8), "negative degrees of freedom")
  data$I1 <- factor(data$I1)
  expect_error(itemrest(data), "Nonnumeric items: I1")
})

test_that("analysis messages are preserved and nonconvergence is flagged", {
  captured <- capture_analysis({ message("important message"); warning("important warning"); 1 })
  expect_equal(captured$messages, "important message")
  expect_equal(captured$warnings, "important warning")
  out <- make_test_efa(names(fixture_data()))
  out$fit_messages <- "rotation did not converge"
  local_mocked_bindings(efa_custom = function(...) out)
  result <- itemrest(fixture_data(), n_factors = 2, reliability = "none", verbose = FALSE)
  expect_equal(result$removal_summary$Convergence, "No")
  expect_equal(result$removal_summary$Messages, "rotation did not converge")
  expect_true(grepl("estimation_message", result$removal_summary$Diagnostics))
  expect_false(result$removal_summary$Candidate)
})

test_that("Howard labels low primary loadings and true cross-loadings separately", {
  out <- list(loadings = matrix(c(.35, .02), 1, dimnames = list("I1", c("F1", "F2"))))
  problems <- identify_problem_items(out, .3, "howard")
  expect_equal(problems$low_loading, "I1")
  expect_length(problems$cross_2, 0L)
  out$loadings[1, ] <- c(.5, .35)
  expect_equal(identify_problem_items(out, .3, "howard")$cross_2, "I1")
})

test_that("congruence matches factor permutations and signs without factorial search", {
  reference <- matrix(c(.8, .7, .6, .05, .1, .05, .05, .1, .05, .7, .8, .6), 6, 2)
  replication <- reference[, c(2, 1)]
  replication[, 1] <- -replication[, 1]
  matched <- matched_congruence(reference, replication)
  expect_equal(matched$order, c(2L, 1L))
  expect_equal(matched$sign, c(1, -1))
  expect_equal(matched$congruence, c(1, 1), ignore_attr = TRUE)
  expect_equal(matched$aligned_loadings, reference)
  expect_equal(matched$min_congruence, 1)
  expect_equal(diag(matched$aligned_matrix), c(1, 1), ignore_attr = TRUE)
  many <- diag(12)
  expect_equal(matched_congruence(many, many[, 12:1])$order, 12:1)
  data <- real_factor_data(n = 300)
  discovery <- itemrest(data[1:150, ], n_factors = 1, verbose = FALSE)
  check <- itemrest_validate(discovery, data[151:300, ])
  expect_true(is.finite(check$validation_summary$Min_Congruence))
})

test_that("fingerprints depend on item names and numeric values, not representations", {
  data <- data.frame(I1 = 1:50, I2 = rep(1:5, 10))
  same <- data
  same$I1 <- same$I1 + 0
  same$I2 <- as.double(same$I2)
  row.names(same) <- paste0("row", seq_len(nrow(same)))
  expect_identical(data_fingerprint(data), data_fingerprint(same))
  same$I1[1] <- 999
  expect_false(identical(data_fingerprint(data), data_fingerprint(same)))
  expect_false(identical(data_fingerprint(data), data_fingerprint(data[, 2:1])))
  data <- real_factor_data(n = 150)
  discovery <- itemrest(data, n_factors = 1, verbose = FALSE)
  same <- data
  same$V1 <- same$V1 + 0
  row.names(same) <- paste0("r", seq_len(nrow(same)))
  expect_silent(boot <- itemrest_bootstrap(discovery, same, n_boot = 2, seed = 7))
  expect_equal(boot$settings$complete_replicates, 2L)
  expect_warning(itemrest_validate(discovery, same), "not independent")
})

test_that("whole-set alpha is labelled and full fits are optional", {
  data <- fixture_data()
  local_mocked_bindings(efa_custom = function(data, ...) make_test_efa(names(data)))
  two <- itemrest(data, n_factors = 2, reliability = "none", verbose = FALSE)
  expect_match(two$removal_summary$Global_Alpha_Interpretation, "descriptive")
  data <- real_factor_data()
  compact <- itemrest(data, n_factors = 1, verbose = FALSE)
  full <- itemrest(data, n_factors = 1, verbose = FALSE, store_fits = TRUE)
  expect_null(compact$initial_efa$efa)
  expect_false(is.null(full$initial_efa$efa))
  expect_equal(compact$removal_summary, full$removal_summary)
})

test_that("all-attempt bootstrap frequencies include unsuccessful searches", {
  data <- fixture_data()
  local_mocked_bindings(efa_custom = function(data, ...) make_test_efa(names(data), low = "I1"))
  discovery <- itemrest(data, n_factors = 2, max_solutions = 1,
    reliability = "none", verbose = FALSE)
  boot <- itemrest_bootstrap(discovery, data, n_boot = 2, seed = 8)
  expect_true(all(is.na(boot$item_stability$Retained_In_Any_Candidate)))
  expect_equal(boot$item_stability$Retained_In_Any_Candidate_All_Attempts, rep(0, 9))
  expect_equal(boot$solution_stability$Candidate_Frequency_All_Attempts, rep(0, nrow(discovery$removal_summary)))
  expect_equal(boot$settings$incomplete_replicates, 2L)
})

test_that("ordinal category thresholds are configurable", {
  data <- data.frame(I1 = rep(1:8, 30), I2 = rep(rep(1:4, each = 2), 30))
  expect_false(is_ordinal_item(data$I1))
  expect_true(is_ordinal_item(data$I1, 8))
  matrix <- cor_matrix_custom(data, "polychoric", "listwise", 8)
  expect_true(all(is.finite(matrix)))
  data <- data.frame(I1 = rep(1:7, 30), I2 = rep(1:5, 42))
  data$I1[1] <- NA
  ordered <- data
  ordered[] <- lapply(ordered, ordered)
  expected <- qgraph::cor_auto(ordered, missing = "pairwise", detectOrdinal = FALSE,
    forcePD = FALSE, verbose = FALSE)
  expect_equal(cor_matrix_custom(data, "polychoric", "pairwise"), expected)
  expect_error(itemrest(data, ordinal_categories = 1), "ordinal_categories")
})

test_that("fa remains the PA default and pc returns component counts", {
  data <- real_factor_data(n = 200, p = 6)
  fa <- itemrest(data, seed = 2, parallel_iterations = 10, reliability = "none", verbose = FALSE)
  pc <- itemrest(data, seed = 2, parallel_iterations = 10, parallel_method = "pc",
    reliability = "none", verbose = FALSE)
  expect_equal(fa$settings$parallel_method, "fa")
  expect_equal(pc$settings$n_factors, pc$parallel_analysis$details$ncomp)
  ordinal <- as.data.frame(lapply(data, function(x) as.integer(cut(x, c(-Inf, -.7, 0, .7, Inf)))))
  original_cor <- cor_matrix_custom
  inputs <- list()
  local_mocked_bindings(cor_matrix_custom = function(data, ...) {
    inputs[[length(inputs) + 1L]] <<- data
    original_cor(data, ...)
  })
  result <- itemrest(ordinal, cor_method = "polychoric", parallel_method = "pc",
    parallel_iterations = 5, seed = 4, reliability = "none", verbose = FALSE)
  expect_equal(length(inputs), 6L) # Once for real observations, once per reference dataset.
  expect_equal(result$parallel_analysis$reference, "independent_column_permutations_same_correlation_method")
  for (input in inputs[-1L]) for (item in names(ordinal))
    expect_equal(table(input[[item]]), table(ordinal[[item]]))
})
