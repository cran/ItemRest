#' Bootstrap the Entire Threshold-Driven Candidate Search
#'
#' @description
#' Resamples observations and reruns discovery for each replicate, holding the
#' original baseline factor count, scoring keys, and screening criteria fixed.
#' Each replicate computes its own correlation matrix once and reuses its
#' submatrices throughout that replicate's removal search.
#' This measures search stability rather than merely refitting a selected set.
#' It does not replace independent validation or estimate probabilities of validity.
#' @param x An itemrest_result defining search settings.
#' @param data Original numeric discovery data, including every discovery item.
#' @param n_boot Number of replicates, default 100.
#' @param seed Integer seed; NULL uses the local date in DDMMYYYY format.
#' @param progress Show a text progress bar; default interactive().
#' @return A list with item_stability, solution_stability, replicates, settings,
#'   and seed. Any/all-candidate frequencies divide by complete successful searches;
#'   searches with no candidates count as zero. Mean candidate retention divides
#'   only by complete replicates with candidates and weights replicates equally.
#'   Failed/incomplete searches are excluded from frequency denominators and
#'   reported separately. Solution frequencies use exact retained sets, not ranks.
#'   Additional *_All_Attempts columns divide by n_boot, conservatively treating
#'   failed and incomplete searches as having no candidate solutions. This is
#'   a sensitivity summary, not an estimate of what those searches would yield.
#' @export
itemrest_bootstrap <- function(x, data, n_boot = 100L, seed = NULL, progress = interactive()) {
  check_itemrest_result(x)
  check_scalar(n_boot, "n_boot", 1, .Machine$integer.max, TRUE)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress))
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  if (!(is.data.frame(data) || is.matrix(data))) stop("data must be a numeric data.frame or matrix.", call. = FALSE)
  data <- as.data.frame(data, check.names = FALSE)
  items <- x$provenance$item_names
  if (any(!items %in% names(data))) stop("Bootstrap data must include all discovery items.", call. = FALSE)
  prepared <- prepare_itemrest_data(data[, items, drop = FALSE], x$settings$missing)
  data <- prepared$data
  if (!identical(data_fingerprint(data), x$provenance$data_fingerprint))
    stop("Bootstrap data must match the original prepared discovery data.", call. = FALSE)
  args <- x$settings[intersect(names(x$settings), names(formals(itemrest)))]
  args$verbose <- FALSE
  args$seed <- NULL
  seed_info <- resolve_itemrest_seed(seed)
  with_itemrest_seed(seed_info$value, {
    any_count <- all_count <- mean_sum <- stats::setNames(numeric(length(items)), items)
    counts <- stats::setNames(numeric(length(x$solution_details)), names(x$solution_details))
    completed <- with_candidates <- 0L
    replications <- vector("list", n_boot)
    # Sampling and analysis seeds are fixed before fitting any replicate.
    bootstrap_rows <- lapply(seq_len(n_boot), function(i) sample.int(nrow(data), nrow(data), replace = TRUE))
    replicate_seeds <- sample.int(.Machine$integer.max, n_boot, replace = TRUE)
    bar <- if (progress) utils::txtProgressBar(min = 0, max = n_boot, style = 3) else NULL
    if (!is.null(bar)) on.exit(close(bar), add = TRUE)
    for (i in seq_len(n_boot)) {
      rows <- bootstrap_rows[[i]]
      replicate_seed <- replicate_seeds[i]
      fitted <- capture_analysis(do.call(itemrest, c(list(data = data[rows, , drop = FALSE], seed = replicate_seed), args)))
      result <- fitted$value
      baseline_failed <- !is.null(result) && result$removal_summary$Analysis_Status[1L] == "skipped"
      failed_branches <- if (is.null(result)) 0L else result$search$failed
      status <- if (!is.null(fitted$error) || baseline_failed || failed_branches > 0L) "failed" else
        if (!result$search$complete) "incomplete" else "complete"
      ids <- if (is.null(result)) character() else result$candidate_solutions$Solution_ID
      if (status == "complete") {
        completed <- completed + 1L
        if (length(ids)) {
          with_candidates <- with_candidates + 1L
          retained <- lapply(ids, function(id) result$solution_details[[id]]$remaining_items)
          fractions <- vapply(items, function(item) mean(vapply(retained, function(set) item %in% set, logical(1))), numeric(1))
          any_count <- any_count + (fractions > 0)
          all_count <- all_count + (fractions == 1)
          mean_sum <- mean_sum + fractions
          for (id in names(x$solution_details)) {
            target <- x$solution_details[[id]]$remaining_items
            counts[id] <- counts[id] + any(vapply(retained, function(set) setequal(set, target), logical(1)))
          }
        }
      }
      error <- if (!is.null(fitted$error)) fitted$error else if (baseline_failed)
        result$removal_summary$Error_Message[1L] else if (failed_branches > 0L)
          result$removal_summary$Error_Message[which(result$removal_summary$Analysis_Status == "failed")[1L]] else ""
      replications[[i]] <- data.frame(Replicate = i, Seed = replicate_seed, Status = status,
        N_Candidates = length(ids), N_Evaluated = if (is.null(result)) 0L else result$search$evaluated,
        N_Failed = failed_branches, N_Skipped = if (is.null(result)) 0L else result$search$skipped,
        Error_Message = error, stringsAsFactors = FALSE)
      if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
    }
    divide <- function(values, denominator) if (denominator) values / denominator else rep(NA_real_, length(values))
    item_stability <- data.frame(Item = items,
      Retained_In_Any_Candidate = unname(divide(any_count, completed)),
      Retained_In_All_Candidates = unname(divide(all_count, completed)),
      Mean_Retention_Among_Candidates = unname(divide(mean_sum, with_candidates)),
      Mean_Removal_Among_Candidates = unname(1 - divide(mean_sum, with_candidates)), stringsAsFactors = FALSE)
    item_stability$Retained_In_Any_Candidate_All_Attempts <- unname(any_count / n_boot)
    item_stability$Retained_In_All_Candidates_All_Attempts <- unname(all_count / n_boot)
    item_stability$Mean_Retention_All_Attempts <- unname(mean_sum / n_boot)
    solution_stability <- data.frame(Solution_ID = names(counts),
      Candidate_Frequency = unname(divide(counts, completed)),
      Candidate_Frequency_All_Attempts = unname(counts / n_boot), stringsAsFactors = FALSE)
    list(item_stability = item_stability, solution_stability = solution_stability,
         replicates = do.call(rbind, replications), seed = seed_info,
         settings = list(n_boot = n_boot, complete_replicates = completed,
                         replicates_with_candidates = with_candidates,
                         failed_replicates = sum(vapply(replications, function(r) r$Status == "failed", logical(1))),
                         incomplete_replicates = sum(vapply(replications, function(r) r$Status == "incomplete", logical(1))),
                         factor_count_policy = "fixed_at_discovery_baseline", n_factors = x$settings$n_factors))
  })
}
