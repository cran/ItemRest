# Correlations are never silently repaired by ItemRest.
#' Calculate a correlation matrix.
#' @keywords internal
cor_matrix_custom <- function(data, method = "polychoric", missing = "pairwise",
                              ordinal_categories = 7L) {
  if (method == "polychoric") {
    data <- as.data.frame(data)
    ordinal <- vapply(data, is_ordinal_item, logical(1), ordinal_categories = ordinal_categories)
    data[ordinal] <- lapply(data[ordinal], ordered)
    invisible(utils::capture.output(result <- qgraph::cor_auto(data,
      missing = if (missing == "pairwise") "pairwise" else "listwise",
      detectOrdinal = FALSE, forcePD = FALSE, verbose = FALSE)))
    result
  } else {
    stats::cor(data, use = if (missing == "pairwise") "pairwise.complete.obs" else "complete.obs",
               method = method)
  }
}

# Each dataset has its own cache. Keep the original, unrepaired matrix so that
# positive-definiteness checks and smoothing remain specific to retained sets.
new_correlation_cache <- function(data, cor_method, missing, ordinal_categories = 7L) {
  cache <- new.env(parent = emptyenv())
  cache$data <- data
  cache$cor_method <- cor_method
  cache$missing <- missing
  cache$ordinal_categories <- ordinal_categories
  cache$counts <- pairwise_counts(data)
  cache$result <- NULL
  cache
}

correlation_counts <- function(data, correlation_cache = NULL) {
  if (is.null(correlation_cache)) return(pairwise_counts(data))
  items <- names(data)
  correlation_cache$counts[items, items, drop = FALSE]
}

cached_correlation_matrix <- function(cache) {
  if (is.null(cache$result)) {
    cache$result <- capture_analysis({
      matrix <- as.matrix(cor_matrix_custom(cache$data, cache$cor_method, cache$missing,
        cache$ordinal_categories))
      p <- ncol(cache$data)
      if (!identical(dim(matrix), c(p, p)))
        stop("Correlation matrix has incorrect dimensions.", call. = FALSE)
      dimnames(matrix) <- list(names(cache$data), names(cache$data))
      matrix
    })
  }
  if (!is.null(cache$result$error)) stop(cache$result$error_condition)
  cache$result$value
}

#' Determine the number of factors using Parallel Analysis.
#' @keywords internal
determine_n_factors <- function(data, cor_method = "pearson", missing = "listwise",
                                pd_action = "fail", parallel_iterations = 100L,
                                correlation_cache = NULL, parallel_method = "fa") {
  correlation <- checked_correlation(data, cor_method, missing, pd_action, correlation_cache)
  counts <- correlation_counts(data, correlation_cache)
  if (cor_method == "polychoric") {
    values <- function(matrix) {
      result <- if (parallel_method == "pc") eigen(matrix, symmetric = TRUE, only.values = TRUE)$values else
        psych::fa(matrix, nfactors = 1L, rotate = "none", fm = "minres",
          warnings = FALSE, smooth = FALSE)$values
      if (any(!is.finite(result))) stop("Parallel-analysis eigenvalues are not finite.", call. = FALSE)
      result
    }
    observed <- values(correlation$matrix)
    ordinal_categories <- if (is.null(correlation_cache)) 7L else correlation_cache$ordinal_categories
    reference <- vapply(seq_len(parallel_iterations), function(i) {
      randomized <- data
      for (item in names(data)) {
        rows <- which(!is.na(data[[item]]))
        randomized[[item]][rows] <- data[[item]][rows[sample.int(length(rows))]]
      }
      cache <- new_correlation_cache(randomized, cor_method, missing, ordinal_categories)
      values(checked_correlation(randomized, cor_method, missing, pd_action, cache)$matrix)
    }, numeric(ncol(data)))
    cutoff <- apply(reference, 1L, stats::quantile, probs = .95)
    first_failure <- which(observed <= cutoff)[1L]
    n_factors <- if (is.na(first_failure)) length(observed) else first_failure - 1L
    return(list(n_factors = n_factors, method = parallel_method,
      reference = "independent_column_permutations_same_correlation_method",
      quantile = .95, observed = observed, cutoff = cutoff, simulated = t(reference)))
  }
  invisible(utils::capture.output(result <- psych::fa.parallel(correlation$matrix, fa = parallel_method, n.iter = parallel_iterations,
                              n.obs = min(counts), plot = FALSE,
                              show.legend = FALSE, main = NULL)))
  list(n_factors = if (parallel_method == "fa") result$nfact else result$ncomp,
       method = parallel_method, reference = "normal_pearson_simulation",
       quantile = .95, details = result)
}

checked_correlation <- function(data, cor_method, missing, pd_action, correlation_cache = NULL) {
  if (any(vapply(data, function(x) length(unique(x[!is.na(x)])) < 2L, logical(1))))
    stop("At least one retained item has no variation.", call. = FALSE)
  counts <- correlation_counts(data, correlation_cache)
  if (min(counts) < 3L) stop("At least three observations are needed for every item pair.", call. = FALSE)
  correlation <- if (is.null(correlation_cache)) {
    as.matrix(cor_matrix_custom(data, cor_method, missing))
  } else {
    items <- names(data)
    cached_correlation_matrix(correlation_cache)[items, items, drop = FALSE]
  }
  if (!identical(dim(correlation), c(ncol(data), ncol(data))) || any(!is.finite(correlation)))
    stop("Correlation matrix contains missing/nonfinite values or has incorrect dimensions.", call. = FALSE)
  correlation <- (correlation + t(correlation)) / 2
  smallest <- min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values)
  positive <- smallest > 1e-8
  adjusted <- FALSE
  if (!positive) {
    if (pd_action == "fail") stop(structure(list(
      message = "Correlation matrix is not positive definite.", call = NULL,
      positive_definite = FALSE, min_eigenvalue = smallest
    ), class = c("itemrest_correlation_error", "error", "condition")))
    correlation <- psych::cor.smooth(correlation)
    adjusted <- TRUE
  }
  list(matrix = correlation, positive_definite = positive,
       min_eigenvalue = smallest, adjusted = adjusted)
}

model_omega <- function(loadings, phi, uniquenesses, index = seq_len(nrow(loadings))) {
  if (length(index) < 2L || any(!is.finite(uniquenesses[index])) ||
      any(uniquenesses[index] < -1e-8)) return(NA_real_)
  common <- loadings[index, , drop = FALSE] %*% phi %*% t(loadings[index, , drop = FALSE])
  signal <- sum(common)
  total <- signal + sum(uniquenesses[index])
  if (!is.finite(total) || total <= 0 || !is.finite(signal) || signal < 0) return(NA_real_)
  value <- signal / total
  if (value < 0 || value > 1 + 1e-8) NA_real_ else min(value, 1)
}

alpha_values <- function(data, correlation) {
  if (ncol(data) < 2L) return(c(raw = NA_real_, standardized = NA_real_, analysis = NA_real_))
  raw <- capture_analysis(psych::alpha(data, check.keys = FALSE, warnings = FALSE,
                                      discrete = FALSE, use = "pairwise"))
  selected <- capture_analysis(psych::alpha(correlation, check.keys = FALSE, warnings = FALSE))
  c(raw = if (is.null(raw$value)) NA_real_ else raw$value$total$raw_alpha,
    standardized = if (is.null(raw$value)) NA_real_ else raw$value$total$std.alpha,
    analysis = if (is.null(selected$value)) NA_real_ else selected$value$total$std.alpha)
}

#' Run an EFA with correlation and numerical diagnostics.
#' @keywords internal
efa_custom <- function(data, n_factors = 1, cor_method = "polychoric", extract = "uls",
                       rotate = "oblimin", missing = "listwise", pd_action = "fail",
                       reliability = "alpha_omega", correlation_cache = NULL) {
  correlation <- checked_correlation(data, cor_method, missing, pd_action, correlation_cache)
  counts <- correlation_counts(data, correlation_cache)
  fitted <- capture_analysis(psych::fa(
    r = correlation$matrix, nfactors = n_factors, rotate = rotate, fm = extract,
    n.obs = min(counts), np.obs = counts, smooth = FALSE
  ))
  if (!is.null(fitted$error)) stop(fitted$error, call. = FALSE)
  efa <- fitted$value
  loadings <- as.matrix(efa$loadings)
  phi <- if (is.null(efa$Phi)) diag(n_factors) else as.matrix(efa$Phi)
  communalities <- diag(loadings %*% phi %*% t(loadings))
  uniquenesses <- 1 - communalities
  alpha <- if (reliability == "none") rep(NA_real_, 3L) else alpha_values(data, correlation$matrix)
  names(alpha) <- c("raw", "standardized", "analysis")
  list(efa = efa, alpha = unname(alpha["raw"]),
       explained_var = mean(communalities), loadings = loadings, phi = phi,
       communalities = communalities, uniquenesses = uniquenesses,
       correlation = correlation, pairwise_n = counts,
       alpha_values = alpha, fit_warnings = fitted$warnings, fit_messages = fitted$messages,
       omega_total = if (reliability == "alpha_omega") model_omega(loadings, phi, uniquenesses) else NA_real_)
}

assess_solution <- function(out, data, n_factors, min_loading, loading_diff,
                            min_items_per_factor, max_factor_correlation, reliability) {
  problems <- identify_problem_items(out, min_loading, loading_diff)
  loadings <- as.matrix(out$loadings)
  phi <- if (is.null(out$phi)) diag(n_factors) else out$phi
  communalities <- diag(loadings %*% phi %*% t(loadings))
  uniquenesses <- 1 - communalities
  finite <- all(is.finite(loadings)) && all(is.finite(phi)) && all(is.finite(communalities))
  notes <- character()
  if (length(problems$low_loading)) notes <- c(notes, "low_loading_items")
  if (length(c(problems$cross_2, problems$cross_3))) notes <- c(notes, "cross_loading_items")
  if (!finite) notes <- c(notes, "nonfinite_estimates")
  heywood <- finite && (any(communalities > 1 + 1e-8) || any(communalities < -1e-8) || any(uniquenesses < -1e-8))
  if (heywood) notes <- c(notes, "inadmissible_variance")
  warnings <- out$fit_warnings
  messages <- out$fit_messages
  message_flags <- estimation_message_flags(messages)
  convergence_warning <- any(grepl("failed.*converg|not.*converg|converg.*(fail|not|limit)|iteration limit|maximum.*iter",
    warnings, ignore.case = TRUE)) || "nonconvergence" %in% message_flags
  reported <- out$efa$converged
  convergence <- if (convergence_warning || identical(reported, FALSE)) "No" else
    if (identical(reported, TRUE)) "Yes" else "Not_reported"
  if (convergence == "No") notes <- c(notes, "nonconvergence")
  max_phi <- if (n_factors > 1L && finite) max(abs(phi[lower.tri(phi)])) else 0
  high_phi <- finite && max_phi >= max_factor_correlation
  if (high_phi) notes <- c(notes, "high_factor_correlation")
  oversized <- finite && any(abs(loadings) > 1 + 1e-8)
  if (oversized) notes <- c(notes, "loading_above_one_requires_review")
  adjusted <- isTRUE(out$correlation$adjusted)
  if (adjusted) notes <- c(notes, "correlation_matrix_smoothed")
  if (length(warnings)) notes <- c(notes, "estimation_warning")
  if (length(message_flags)) notes <- c(notes, "estimation_message", message_flags)
  assignment <- if (finite) max.col(abs(loadings), ties.method = "first") else rep(NA_integer_, nrow(loadings))
  names(assignment) <- rownames(loadings)
  qualified <- assignment
  qualified[names(qualified) %in% all_flagged_items(problems)] <- NA_integer_
  factor_counts <- tabulate(qualified, nbins = n_factors)
  supported <- all(factor_counts >= min_items_per_factor)
  if (!supported) notes <- c(notes, "too_few_qualifying_items_per_factor")
  dof <- ncol(data) * (ncol(data) - 1) / 2 - ncol(data) * n_factors + n_factors * (n_factors - 1) / 2
  phi_positive <- finite && all(eigen(phi, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  if (!phi_positive) notes <- c(notes, "invalid_factor_correlation_matrix")
  admissible <- finite && !heywood && convergence != "No" && dof >= 0 &&
    (is.null(out$correlation) || isTRUE(out$correlation$positive_definite) || adjusted) &&
    phi_positive
  review <- high_phi || oversized || adjusted || length(warnings) > 0L || length(message_flags) > 0L
  candidate <- admissible && !review && supported && length(all_flagged_items(problems)) == 0L
  factor_reliability <- lapply(seq_len(n_factors), function(factor) {
    index <- which(qualified == factor)
    values <- if (length(index) >= 2L && reliability != "none") {
      alpha_values(data[, index, drop = FALSE], out$correlation$matrix[index, index, drop = FALSE])
    } else rep(NA_real_, 3L)
    data.frame(Factor = colnames(loadings)[factor], N_Items = length(index),
               Items = item_label(names(qualified)[index]),
               Raw_Alpha = unname(values[1L]), Standardized_Alpha = unname(values[2L]),
               Analysis_Correlation_Alpha = unname(values[3L]),
               Omega_Total = if (reliability == "alpha_omega" && admissible)
                 model_omega(loadings, phi, uniquenesses, index) else NA_real_,
               stringsAsFactors = FALSE)
  })
  list(problems = problems, assignment = assignment, qualified_assignment = qualified,
       factor_counts = factor_counts, factor_reliability = do.call(rbind, factor_reliability),
       candidate = candidate, admissible = admissible, review = review,
       convergence = convergence, max_phi = max_phi, dof = dof,
       communalities = communalities, uniquenesses = uniquenesses, notes = unique(notes))
}
