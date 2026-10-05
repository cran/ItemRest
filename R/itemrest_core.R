#' Check the Howard loading criteria.
#' @keywords internal
howard <- function(primary, secondary) {
  tolerance <- 1e-12
  !(primary >= 0.40 - tolerance && secondary < 0.30 - tolerance &&
      primary - secondary >= 0.20 - tolerance)
}

#' Calculate basic descriptive statistics.
#' @keywords internal
descriptive_stats <- function(data) {
  list(n_items = ncol(data), n_obs = nrow(data), min_value = min(data, na.rm = TRUE),
       max_value = max(data, na.rm = TRUE))
}

#' Identify low-loading and cross-loading items.
#' @keywords internal
identify_problem_items <- function(efa_res, min_loading = 0.30, loading_diff = 0.10) {
  loadings <- abs(as.matrix(efa_res$loadings))
  result <- empty_problem_items()
  if (!nrow(loadings) || !ncol(loadings)) return(result)
  for (i in seq_len(nrow(loadings))) {
    values <- loadings[i, ]
    item <- rownames(loadings)[i]
    primary_min <- if (identical(loading_diff, "howard")) max(min_loading, .40) else min_loading
    if (any(!is.finite(values)) || max(values) < primary_min - 1e-12) {
      result$low_loading <- c(result$low_loading, item)
      next
    }
    if (length(values) > 1L) {
      sorted <- sort(values, decreasing = TRUE)
      cross <- if (identical(loading_diff, "howard")) howard(sorted[1L], sorted[2L]) else
        sorted[1L] - sorted[2L] < loading_diff - 1e-12
      if (cross) {
        category <- if (sum(values >= min_loading) >= 3L) "cross_3" else "cross_2"
        result[[category]] <- c(result[[category]], item)
      }
    }
  }
  result
}

evaluate_item_set <- function(data, n_factors, cor_method, extract, rotate,
                              min_loading, loading_diff, options, correlation_cache = NULL) {
  p <- ncol(data)
  dof <- p * (p - 1) / 2 - p * n_factors + n_factors * (n_factors - 1) / 2
  if (p <= n_factors || dof < 0) return(list(
    status = "skipped", error = "Too few items or negative EFA degrees of freedom.",
    warnings = character(), messages = character(), out = NULL, assessment = NULL
  ))
  analysis <- capture_analysis(do.call(efa_custom, c(list(
    data = data, n_factors = n_factors, cor_method = cor_method, extract = extract,
    rotate = rotate, correlation_cache = correlation_cache
  ), options[c("missing", "pd_action", "reliability")])))
  if (!is.null(analysis$error)) return(list(
    status = "failed", error = analysis$error, warnings = analysis$warnings, messages = analysis$messages,
    out = NULL, assessment = NULL, error_condition = analysis$error_condition
  ))
  assessed <- capture_analysis(do.call(assess_solution, c(list(
    out = analysis$value, data = data, n_factors = n_factors,
    min_loading = min_loading, loading_diff = loading_diff
  ), options[c("min_items_per_factor", "max_factor_correlation", "reliability")])))
  if (!is.null(assessed$error)) return(list(
    status = "failed", error = assessed$error, warnings = c(analysis$warnings, assessed$warnings),
    messages = unique(c(analysis$messages, assessed$messages)),
    out = analysis$value, assessment = NULL
  ))
  list(status = if (assessed$value$admissible) "ok" else "inadmissible",
       error = NULL, warnings = unique(c(analysis$warnings, analysis$value$fit_warnings, assessed$warnings)),
       messages = unique(c(analysis$messages, analysis$value$fit_messages, assessed$messages)),
       out = analysis$value, assessment = assessed$value)
}

solution_row <- function(id, state, base_items, evaluated, data, correlation_cache = NULL) {
  remaining <- setdiff(base_items, state$removed)
  out <- evaluated$out
  assessment <- evaluated$assessment
  assessed <- !is.null(assessment)
  loadings <- if (assessed) out$loadings else matrix(numeric(), 0, 0)
  valid <- loadings[is.finite(loadings) & abs(loadings) >= state$min_loading]
  counts <- correlation_counts(data, correlation_cache)
  alpha <- out$alpha_values
  data.frame(
    Solution_ID = id, Baseline = state$iteration == 0L, Iteration = state$iteration,
    Removed_This_Step = item_label(state$step), Removed_Items = item_label(state$removed),
    N_Removed = length(state$removed), Remaining_Items = item_label(remaining),
    N_Remaining = length(remaining),
    Total_Explained_Var = if (assessed) out$explained_var else NA_real_,
    Factor_Loading_Range = if (length(valid)) paste0(formatC(min(valid), digits = 2, format = "f"),
      "/", formatC(max(valid), digits = 2, format = "f")) else "N/A",
    Absolute_Loading_Range = if (length(valid)) paste0(formatC(min(abs(valid)), digits = 2, format = "f"),
      "/", formatC(max(abs(valid)), digits = 2, format = "f")) else "N/A",
    Cronbachs_Alpha = if (!is.null(out$alpha)) out$alpha else NA_real_,
    Global_Alpha_Interpretation = if (!assessed) "Not_assessed" else if (ncol(loadings) > 1L)
      "whole_set_descriptive_not_unidimensional_reliability" else "whole_set",
    Standardized_Alpha = if (!is.null(alpha)) unname(alpha["standardized"]) else NA_real_,
    Analysis_Correlation_Alpha = if (!is.null(alpha)) unname(alpha["analysis"]) else NA_real_,
    Omega_Total = if (assessed && assessment$admissible && !is.null(out$omega_total)) out$omega_total else NA_real_,
    Cross_Loading = if (!assessed) "Not_assessed" else if (length(c(assessment$problems$cross_2,
      assessment$problems$cross_3))) "Yes" else "No",
    Low_Loading = if (!assessed) "Not_assessed" else if (length(assessment$problems$low_loading)) "Yes" else "No",
    Problem_Items = if (assessed) item_label(all_flagged_items(assessment$problems)) else "Not_assessed",
    Analysis_Status = evaluated$status,
    Candidate = assessed && assessment$candidate,
    Admissible = if (assessed) assessment$admissible else NA,
    Review_Required = if (assessed) assessment$review else NA,
    N_Obs = if (length(counts)) min(counts) else NA_real_,
    N_Obs_Max = if (length(counts)) max(counts) else NA_real_, N_Rows = nrow(data),
    Positive_Definite = if (!is.null(out$correlation)) out$correlation$positive_definite else
      if (!is.null(evaluated$error_condition$positive_definite)) evaluated$error_condition$positive_definite else NA,
    Min_Eigenvalue = if (!is.null(out$correlation)) out$correlation$min_eigenvalue else
      if (!is.null(evaluated$error_condition$min_eigenvalue)) evaluated$error_condition$min_eigenvalue else NA_real_,
    Correlation_Adjusted = if (!is.null(out$correlation)) out$correlation$adjusted else NA,
    Convergence = if (assessed) assessment$convergence else "Not_assessed",
    Max_Communality = if (assessed) max(assessment$communalities) else NA_real_,
    Min_Uniqueness = if (assessed) min(assessment$uniquenesses) else NA_real_,
    Max_Abs_Loading = if (assessed) max(abs(loadings)) else NA_real_,
    Max_Factor_Correlation = if (assessed) assessment$max_phi else NA_real_,
    Min_Factor_Items = if (assessed) min(assessment$factor_counts) else NA_integer_,
    Degrees_Of_Freedom = if (assessed) assessment$dof else NA_real_,
    Diagnostics = if (assessed) paste(assessment$notes, collapse = "; ") else evaluated$error,
    Warnings = paste(evaluated$warnings, collapse = "; "),
    Messages = paste(evaluated$messages, collapse = "; "),
    Error_Message = if (is.null(evaluated$error)) "" else evaluated$error,
    stringsAsFactors = FALSE
  )
}

rank_solution_table <- function(table, rank_by) {
  if (rank_by == "n_removed") table <- table[order(!table$Baseline, table$N_Removed), , drop = FALSE]
  if (rank_by == "explained_variance")
    table <- table[order(!table$Baseline, -table$Total_Explained_Var, na.last = TRUE), , drop = FALSE]
  row.names(table) <- NULL
  table
}

#' Iteratively evaluate removal strategies and retain failed attempts.
#' @keywords internal
test_removals <- function(data, base_items, n_factors, cor_method, extract, rotate,
                          min_loading, loading_diff,
                          options = analysis_options("listwise", 3L, 0.85, "fail", "alpha_omega"),
                          retain_items = character(), max_solutions = 10000L,
                          correlation_cache = NULL, store_fits = FALSE) {
  if (is.null(correlation_cache))
    correlation_cache <- new_correlation_cache(data, cor_method, options$missing)
  queue <- list(list(removed = character(), step = character(), iteration = 0L, parent = NA_character_))
  scheduled <- new.env(hash = TRUE, parent = emptyenv())
  assign("items:", TRUE, envir = scheduled)
  rows <- details <- list()
  truncated <- FALSE
  enqueue <- function(removed, step, iteration, parent) {
    removed <- base_items[base_items %in% removed]
    if (any(removed %in% retain_items)) return(invisible(NULL))
    key <- paste0("items:", paste(which(base_items %in% removed), collapse = ","))
    if (exists(key, envir = scheduled, inherits = FALSE)) return(invisible(NULL))
    if (length(queue) >= max_solutions) {
      truncated <<- TRUE
      return(invisible(NULL))
    }
    assign(key, TRUE, envir = scheduled)
    queue[[length(queue) + 1L]] <<- list(removed = removed, step = step,
                                      iteration = iteration, parent = parent)
    invisible(NULL)
  }
  enqueue_combinations <- function(items, state, parent) {
    items <- setdiff(items, retain_items)
    for (size in seq_along(items)) {
      index <- seq_len(size)
      repeat {
        step <- items[index]
        enqueue(c(state$removed, step), step, state$iteration + 1L, parent)
        if (truncated) return(invisible(NULL))
        positions <- which(index < length(items) - size + seq_len(size))
        if (!length(positions)) break
        j <- max(positions)
        index[j] <- index[j] + 1L
        if (j < size) index[(j + 1L):size] <- index[j] + seq_len(size - j)
      }
    }
    invisible(NULL)
  }
  i <- 1L
  while (i <= length(queue)) {
    state <- queue[[i]]
    state$min_loading <- min_loading
    id <- sprintf("S%05d", i)
    remaining <- setdiff(base_items, state$removed)
    subset <- data[, remaining, drop = FALSE]
    evaluated <- evaluate_item_set(subset, n_factors, cor_method, extract, rotate,
                                   min_loading, loading_diff, options, correlation_cache)
    rows[[i]] <- solution_row(id, state, base_items, evaluated, subset, correlation_cache)
    if (!store_fits && !is.null(evaluated$out)) evaluated$out$efa <- NULL
    details[[id]] <- list(remaining_items = remaining, removed_items = state$removed,
                          removed_this_step = state$step, parent_id = state$parent,
                          analysis = evaluated$out, assessment = evaluated$assessment,
                          status = evaluated$status, error = evaluated$error)
    if (!is.null(evaluated$assessment))
      enqueue_combinations(all_flagged_items(evaluated$assessment$problems), state, id)
    i <- i + 1L
  }
  list(summary = do.call(rbind, rows), details = details, truncated = truncated,
       evaluated = length(rows), scheduled = length(queue))
}

#' Identify Candidate Item Sets for Exploratory Factor Analysis
#'
#' @description
#' Reassesses remaining items after every removal and evaluates combinations of
#' flagged items. This threshold-driven search in one sample does not establish
#' an optimal, valid, or independently replicated measurement solution.
#' Factor count is determined once at baseline and stays fixed.
#' The original correlation matrix and observation counts are computed once
#' for all baseline items, after missing-data handling and scoring keys. Each
#' retained set uses the corresponding rows and columns. Positive-definiteness
#' checks and any requested smoothing are applied separately to each subset.
#' @param data Numeric data.frame or matrix, with unique item names.
#' @param cor_method "pearson", "spearman", "kendall", or "polychoric". The last
#'   uses qgraph mixed correlations after detecting integer-valued ordinal items
#'   up to ordinal_categories observed categories (excluding missing values).
#' @param n_factors Fixed factor count, or NULL for baseline parallel analysis.
#' @param extract Extraction method passed to psych::fa; default "uls".
#' @param rotate Rotation passed to psych::fa; default "oblimin".
#' @param min_loading Minimum absolute primary loading; default 0.30.
#' @param loading_diff Minimum primary-secondary difference, or "howard".
#' @param missing "listwise" (default), "pairwise", or "fail". Listwise handling
#'   is applied once across all baseline items, fixing observations across sets.
#'   Pairwise counts are recorded; their minimum is used as conservative n.obs.
#'   This does not resolve missing-data bias.
#' @param min_items_per_factor Minimum qualifying primary items per factor;
#'   default 3. Flagged items do not count toward this threshold.
#' @param max_factor_correlation Absolute factor correlation requiring review;
#'   default 0.85, a configurable screening threshold rather than a validity rule.
#' @param pd_action "fail" or "smooth" for nonpositive definite correlations.
#'   Smoothed sets require review and are excluded from candidate screening.
#' @param reliability "alpha_omega", "alpha", or "none". Raw, standardized,
#'   and selected-correlation alpha are distinguished. Model omega total uses
#'   sum(L Phi L') / (sum(L Phi L') + sum(uniquenesses)) for a unit-weighted
#'   standardized sum. It includes all common factors, is not omega hierarchical,
#'   and does not establish unidimensionality. No items are automatically reversed.
#' @param rank_by "none" (discovery order), "n_removed" (ascending), or
#'   "explained_variance" (descending). No winner is selected. Explained variance
#'   is mean model communality within the retained set; values across different
#'   sets do not establish superiority.
#' @param retain_items Item names protected from removal for content reasons.
#'   Protection does not waive screening criteria. A flagged protected item can
#'   prevent all candidate solutions; inspect Problem_Items and content_decisions.
#' @param item_reasons Named character vector documenting item decisions.
#' @param max_solutions Maximum queued/evaluated sets including baseline;
#'   default 10000. A bounded search is explicitly labelled incomplete.
#' @param seed Integer seed. NULL uses the execution host's local calendar date
#'   in DDMMYYYY format (e.g. 05102026 is used numerically as 5102026). Both the
#'   formatted label and numeric seed are recorded. RNG state is restored on exit.
#' @param keys Named numeric vector of 1 or -1 for specified items. Minus one
#'   explicitly reverses scoring by negating values. No automatic reverse scoring
#'   occurs; offsets do not affect correlations or reliability.
#' @param parallel_iterations Number of parallel-analysis replications; default 100.
#' @param parallel_method "fa" (default) for reduced-matrix eigenvalues based on
#'   a one-factor minres fit, or "pc" for full-matrix component eigenvalues.
#'   Polychoric analysis permutes observed values within each item, preserving
#'   categories and missing positions, and computes the same correlation type
#'   for each reference sample. The 95th percentile is used for comparison.
#' @param ordinal_categories Maximum number of integer-valued categories detected
#'   as ordinal by the polychoric backend; default 7.
#' @param store_fits Keep full psych EFA fit objects in solution_details; default
#'   FALSE retains loadings, Phi, communalities, reliability, and diagnostics.
#' @param verbose Print a brief search summary; default TRUE.
#' @return An itemrest_result containing candidate_solutions, removal_summary
#'   (baseline, failures, and skipped sets), solution_details (EFA, diagnostics,
#'   assignments, factor reliability, and first-discovered paths), initial_efa,
#'   problem_items, descriptive_stats, settings, search, content_decisions,
#'   data_summary, provenance, correlation_matrix (the original unsmoothed
#'   baseline matrix), and pairwise_n. correlation_matrix is NULL when no
#'   correlation could be computed. candidate_solutions may have zero rows.
#'   correlation_warnings and correlation_messages belong to the source matrix;
#'   they are recorded once and not attributed to every retained set. Estimator
#'   messages are recorded per set, but informational messages do not require
#'   review; explicit nonconvergence, unavailable rotation, variance problems,
#'   and matrix repair reports do.
#'   parallel_analysis records automatic factor-count diagnostics. Full EFA fit
#'   objects are included only when store_fits = TRUE.
#'   Screening eligibility is not substantive validity. The old optimal_strategy
#'   field is replaced by candidate_solutions.
#'   search$complete is FALSE for a limit or numerical branch failure; skipped
#'   underidentified sets remain documented terminal attempts.
#' @export
#' @examples
#' set.seed(4)
#' f <- rnorm(150)
#' d <- data.frame(I1 = f + rnorm(150), I2 = f + rnorm(150),
#'                 I3 = f + rnorm(150), I4 = f + rnorm(150))
#' result <- itemrest(d, n_factors = 1, seed = 10, verbose = FALSE)
#' print(result)
itemrest <- function(data, cor_method = "pearson", n_factors = NULL, extract = "uls",
                     rotate = "oblimin", min_loading = 0.30, loading_diff = 0.10,
                     missing = c("listwise", "pairwise", "fail"), min_items_per_factor = 3L,
                     max_factor_correlation = 0.85, pd_action = c("fail", "smooth"),
                     reliability = c("alpha_omega", "alpha", "none"),
                     rank_by = c("none", "n_removed", "explained_variance"),
                     retain_items = character(), item_reasons = NULL, max_solutions = 10000L,
                     seed = NULL, verbose = TRUE, keys = NULL, parallel_iterations = 100L,
                     ordinal_categories = 7L, store_fits = FALSE,
                     parallel_method = c("fa", "pc")) {
  missing <- match.arg(missing)
  pd_action <- match.arg(pd_action)
  reliability <- match.arg(reliability)
  rank_by <- match.arg(rank_by)
  parallel_method <- match.arg(parallel_method)
  cor_method <- match.arg(tolower(cor_method), c("pearson", "spearman", "kendall", "polychoric"))
  methods <- validate_efa_methods(extract, rotate)
  extract <- methods$extract
  rotate <- methods$rotate
  check_scalar(ordinal_categories, "ordinal_categories", 2L, .Machine$integer.max, TRUE)
  if (!is.logical(store_fits) || length(store_fits) != 1L || is.na(store_fits))
    stop("store_fits must be TRUE or FALSE.", call. = FALSE)
  check_scalar(min_loading, "min_loading", 0, 1)
  if (!identical(loading_diff, "howard")) check_scalar(loading_diff, "loading_diff", 0, 1)
  check_scalar(min_items_per_factor, "min_items_per_factor", 1, .Machine$integer.max, TRUE)
  check_scalar(max_factor_correlation, "max_factor_correlation", .Machine$double.eps, 1)
  check_scalar(max_solutions, "max_solutions", 1, .Machine$integer.max, TRUE)
  check_scalar(parallel_iterations, "parallel_iterations", 1, .Machine$integer.max, TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  prepared <- prepare_itemrest_data(data, missing)
  data <- prepared$data
  fingerprint <- data_fingerprint(data)
  keyed <- apply_item_keys(data, keys)
  data <- keyed$data
  if (!is.character(retain_items) || anyNA(retain_items) || any(!retain_items %in% names(data)))
    stop("retain_items must contain known item names.", call. = FALSE)
  if (!is.null(item_reasons) && (!is.character(item_reasons) || is.null(names(item_reasons)) ||
      anyNA(item_reasons) || anyDuplicated(names(item_reasons)) || any(!names(item_reasons) %in% names(data))))
    stop("item_reasons must be a named character vector for known items.", call. = FALSE)
  auto <- is.null(n_factors)
  seed_info <- resolve_itemrest_seed(seed)
  seed <- seed_info$value
  with_itemrest_seed(seed, {
    correlation_cache <- new_correlation_cache(data, cor_method, missing, ordinal_categories)
    factor_determination <- if (auto) capture_analysis(determine_n_factors(data, cor_method,
      missing, pd_action, parallel_iterations, correlation_cache, parallel_method)) else NULL
    if (auto) {
      if (!is.null(factor_determination$error)) stop("Baseline parallel analysis failed: ", factor_determination$error, call. = FALSE)
      n_factors <- factor_determination$value$n_factors
      if (is.na(n_factors) || n_factors < 1L)
        stop("Parallel analysis did not suggest a positive factor count; inspect the data or supply n_factors.", call. = FALSE)
    }
    check_scalar(n_factors, "n_factors", 1, ncol(data) - 1L, TRUE)
    p <- ncol(data)
    if (p * (p - 1) / 2 - p * n_factors + n_factors * (n_factors - 1) / 2 < 0)
      stop("Baseline EFA has negative degrees of freedom; reduce n_factors or supply more items.", call. = FALSE)
    options <- analysis_options(missing, min_items_per_factor, max_factor_correlation, pd_action, reliability)
    search <- test_removals(data, names(data), n_factors, cor_method, extract, rotate,
                            min_loading, loading_diff, options, retain_items, max_solutions,
                            correlation_cache, store_fits)
    summary <- rank_solution_table(search$summary, rank_by)
    candidates <- summary[summary$Candidate, , drop = FALSE]
    row.names(candidates) <- NULL
    baseline <- search$details[["S00001"]]
    n_failed <- sum(summary$Analysis_Status == "failed")
    n_skipped <- sum(summary$Analysis_Status == "skipped")
    problems <- if (is.null(baseline$assessment)) empty_problem_items() else baseline$assessment$problems
    settings <- c(list(n_factors = n_factors, cor_method = cor_method, extract = extract,
      rotate = rotate, min_loading = min_loading, loading_diff = loading_diff), options,
      list(rank_by = rank_by, retain_items = unique(retain_items), item_reasons = item_reasons,
           keys = keyed$keys, parallel_iterations = parallel_iterations,
           parallel_method = parallel_method,
           ordinal_categories = ordinal_categories, store_fits = store_fits,
           max_solutions = max_solutions, seed = seed, seed_label = seed_info$label,
           seed_source = seed_info$source, auto_n_factors = auto,
           factor_count_policy = "fixed_at_baseline",
           correlation_policy = "baseline_submatrix", verbose = verbose))
    reasons <- rep(NA_character_, ncol(data))
    names(reasons) <- names(data)
    reasons[names(item_reasons)] <- item_reasons
    result <- structure(list(
      descriptive_stats = descriptive_stats(data), initial_efa = baseline$analysis,
      problem_items = problems, all_problem_items_combined = all_flagged_items(problems),
      candidate_solutions = candidates, removal_summary = summary,
      solution_details = search$details, settings = settings,
      search = list(evaluated = search$evaluated, complete = !search$truncated && n_failed == 0L,
                    limit_reached = search$truncated, failed = n_failed, skipped = n_skipped,
                    termination = if (search$truncated) "max_solutions" else if (n_failed) "numerical_failures" else "completed",
                    scope = "threshold_reachable_item_sets", max_solutions = max_solutions),
      content_decisions = data.frame(Item = names(data), Protected = names(data) %in% retain_items,
                                     Reason = unname(reasons), stringsAsFactors = FALSE),
      data_summary = prepared$summary,
      correlation_matrix = correlation_cache$result$value,
      correlation_warnings = as.character(correlation_cache$result$warnings),
      correlation_messages = as.character(correlation_cache$result$messages),
      pairwise_n = correlation_cache$counts,
      parallel_analysis = factor_determination$value,
      provenance = c(package_provenance(), list(seed = seed, seed_label = seed_info$label,
                      seed_source = seed_info$source, item_names = names(data),
                      ordinal_items = if (cor_method == "polychoric")
                        names(data)[vapply(data, is_ordinal_item, logical(1),
                          ordinal_categories = ordinal_categories)] else character(),
                      data_fingerprint = fingerprint,
                      validation_fingerprint = data_fingerprint(prepared$data, ignore_row_order = TRUE),
                      parallel_analysis_warnings = factor_determination$warnings))
    ), class = "itemrest_result")
    if (verbose) {
      message("Evaluated ", nrow(summary), " item sets; ", nrow(candidates), " candidate solutions.")
      if (search$truncated) message("Search incomplete: max_solutions reached.")
    }
    result
  })
}

#' Print Candidate Solutions from ItemRest
#' @param x An itemrest_result.
#' @param report "candidates" (default) or "all". The old "optimal" value is a
#'   deprecated alias for "candidates".
#' @param ... Unused.
#' @return x, invisibly.
#' @export
print.itemrest_result <- function(x, report = c("candidates", "all"), ...) {
  if (identical(report, "optimal")) {
    warning("report = 'optimal' is deprecated; use 'candidates'.", call. = FALSE)
    report <- "candidates"
  }
  report <- match.arg(report)
  cat("\nItemRest: Candidate Solutions\n")
  cat("Factor count:", x$settings$n_factors, "(fixed at baseline); ordering:", x$settings$rank_by, "\n")
  cat("Seed:", x$settings$seed_label, "\n")
  cat("Search:", if (x$search$complete) "complete within threshold-driven scope" else
    paste0("incomplete (", x$search$termination, ")"), "\n")
  table <- if (report == "all") x$removal_summary else x$candidate_solutions
  if (report == "candidates" && !"S00001" %in% table$Solution_ID) {
    cat("No-removal baseline:\n")
    baseline <- x$removal_summary[x$removal_summary$Baseline, , drop = FALSE]
    print(baseline[, c("Solution_ID", "Total_Explained_Var", "Analysis_Status", "Candidate", "Diagnostics")], row.names = FALSE)
    cat("Candidate solutions:\n")
  }
  if (!nrow(table)) cat("No candidate solutions meet the screening criteria. Inspect the baseline and diagnostics with report = 'all'.\n")
  else {
    columns <- c("Solution_ID", "Baseline", "Removed_Items", "N_Remaining", "Total_Explained_Var",
                 "Analysis_Status", "Candidate", "Diagnostics")
    print(table[, columns, drop = FALSE], row.names = FALSE)
  }
  cat("Candidates require content review and independent evaluation. Explained variance across different item sets is descriptive.\n")
  cat("Let algorithms be your compass, not your captain.\n")
  invisible(x)
}
