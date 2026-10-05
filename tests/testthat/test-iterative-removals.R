test_that("new cross-loading and low-loading items are found in later rounds", {
  data <- fixture_data()
  local_mocked_bindings(efa_custom = function(data, ...) {
    remaining <- names(data)
    if ("I1" %in% remaining) return(make_test_efa(remaining, low = "I1"))
    if ("I2" %in% remaining) return(make_test_efa(remaining, cross = "I2"))
    if ("I3" %in% remaining) return(make_test_efa(remaining, low = "I3"))
    make_test_efa(remaining)
  })
  result <- itemrest(data, n_factors = 2, reliability = "none", verbose = FALSE)
  path <- result$removal_summary
  expect_equal(result$all_problem_items_combined, "I1")
  expect_equal(path$Iteration, 0:3)
  expect_equal(path$Problem_Items, c("I1", "I2", "I3", "None"))
  expect_equal(path$Removed_This_Step, c("None", "I1", "I2", "I3"))
  expect_equal(path$N_Removed, 0:3)
  expect_equal(path$N_Remaining, 9:6)
  expect_equal(result$candidate_solutions$Removed_Items, "I1-I2-I3")
  expect_false("optimal_strategy" %in% names(result))
  expect_match(paste(capture.output(print(result)), collapse = "\n"), "No-removal baseline")
  expect_equal(result$solution_details[["S00004"]]$parent_id, "S00003")
})

test_that("converging paths evaluate each remaining item set once", {
  data <- fixture_data()
  calls <- character()
  local_mocked_bindings(efa_custom = function(data, ...) {
    remaining <- names(data)
    calls <<- c(calls, paste(remaining, collapse = ","))
    low <- intersect(if (length(remaining) == 9L) c("I1", "I2") else c("I1", "I2", "I3"), remaining)
    make_test_efa(remaining, low = low)
  })
  result <- itemrest(data, n_factors = 2, reliability = "none", verbose = FALSE)
  expect_equal(length(calls), length(unique(calls)))
  expect_equal(nrow(result$removal_summary), 7L)
  expect_equal(sum(result$removal_summary$Removed_Items == "I1-I2-I3"), 1L)
})

test_that("failed branches are recorded without aborting other branches", {
  data <- fixture_data()
  local_mocked_bindings(efa_custom = function(data, ...) {
    removed <- setdiff(names(fixture_data()), names(data))
    if (identical(removed, "I1")) stop("deliberate branch failure")
    make_test_efa(names(data), low = if (!length(removed)) c("I1", "I2") else character())
  })
  result <- itemrest(data, n_factors = 2, reliability = "none", verbose = FALSE)
  expect_equal(sum(result$removal_summary$Analysis_Status == "failed"), 1L)
  expect_true("I2" %in% result$candidate_solutions$Removed_Items)
  expect_match(result$removal_summary$Error_Message[2L], "deliberate branch failure")
  expect_equal(result$removal_summary$Cross_Loading[2L], "Not_assessed")
})

test_that("infeasible states are recorded and never fitted", {
  local_mocked_bindings(efa_custom = function(...) stop("must not be called"))
  expect_error(itemrest(fixture_data(p = 3), n_factors = 2, verbose = FALSE), "negative degrees of freedom")
  data <- fixture_data(p = 3)
  result <- test_removals(data, names(data), 2, "pearson", "uls", "oblimin", .3, .1)
  expect_equal(result$summary$Analysis_Status, "skipped")
})

test_that("all screened candidates are retained and ranking is explicit", {
  local_mocked_bindings(efa_custom = function(data, ...) {
    removed <- setdiff(names(fixture_data()), names(data))
    out <- make_test_efa(names(data), low = if (!length(removed)) c("I1", "I2") else character())
    out$explained_var <- if (!length(removed)) 0.99 else if (identical(removed, "I2")) 0.8 else 0.6
    out
  })
  default <- itemrest(fixture_data(), n_factors = 2, reliability = "none", verbose = FALSE)
  ranked <- itemrest(fixture_data(), n_factors = 2, reliability = "none", verbose = FALSE,
                    rank_by = "explained_variance")
  expect_equal(nrow(default$candidate_solutions), 3L)
  expect_equal(default$candidate_solutions$Removed_Items[1L], "I1")
  expect_equal(ranked$candidate_solutions$Removed_Items[1L], "I2")
  expect_true(ranked$removal_summary$Baseline[1L])
  expect_warning(print(default, report = "optimal"), "deprecated")
  expect_error(print(default, report = "typo"), "arg")
})

test_that("baseline is always available and counts as a candidate when eligible", {
  local_mocked_bindings(efa_custom = function(data, ...) make_test_efa(names(data)))
  result <- itemrest(fixture_data(), n_factors = 2, reliability = "none", verbose = FALSE)
  expect_equal(nrow(result$removal_summary), 1L)
  expect_equal(nrow(result$candidate_solutions), 1L)
  expect_true(result$removal_summary$Baseline)
})

test_that("search limits and protected items are respected", {
  local_mocked_bindings(efa_custom = function(data, ...) {
    make_test_efa(names(data), low = intersect(names(data), c("I1", "I2", "I3")))
  })
  limited <- itemrest(fixture_data(), n_factors = 2, reliability = "none", verbose = FALSE, max_solutions = 2)
  protected <- itemrest(fixture_data(), n_factors = 2, reliability = "none", verbose = FALSE,
                        retain_items = "I1", item_reasons = c(I1 = "Required content coverage"))
  expect_false(limited$search$complete)
  expect_equal(nrow(limited$removal_summary), 2L)
  expect_true(all(vapply(protected$solution_details, function(s) "I1" %in% s$remaining_items, logical(1))))
  expect_equal(nrow(protected$candidate_solutions), 0L)
  expect_equal(protected$content_decisions$Reason[1L], "Required content coverage")
})

test_that("Howard failures below the old secondary threshold are captured", {
  out <- list(loadings = matrix(c(0.35, 0.10), nrow = 1, dimnames = list("I1", c("F1", "F2"))))
  expect_equal(identify_problem_items(out, 0.30, "howard")$low_loading, "I1")
  expect_length(identify_problem_items(out, 0.30, "howard")$cross_2, 0L)
  expect_equal(identify_problem_items(out, 0.30, 0.10)$cross_2, character())
  expect_equal(identify_problem_items(out, 0.40, 0.10)$low_loading, "I1")
})

test_that("decimal threshold boundaries are not changed by rounding error", {
  out <- list(loadings = matrix(c(.5, .4), nrow = 1, dimnames = list("I1", c("F1", "F2"))))
  expect_length(identify_problem_items(out, .3, .1)$cross_2, 0L)
  expect_false(howard(.45, .25))
})
