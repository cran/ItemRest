# ===================================================================
# NOTE ON PACKAGE DEPENDENCIES
# ===================================================================
# Required packages (e.g., psych, gtools, utils, stats) should be
# listed under the 'Imports:' field in the package's DESCRIPTION file.
# ===================================================================


# ===================================================================
# HELPER (INTERNAL) FUNCTIONS
# ===================================================================

#' Calculate basic descriptive statistics for a dataset.
#' @keywords internal
descriptive_stats <- function(data) {
  list(
    n_items = ncol(data),
    n_obs = nrow(data),
    min_value = min(data, na.rm = TRUE),
    max_value = max(data, na.rm = TRUE)
  )
}

#' Calculate a correlation matrix.
#' @keywords internal
cor_matrix_custom <- function(data, method = "polychoric") {
  if (method == "polychoric") {
    cor_mat <- suppressMessages(qgraph::cor_auto(data))
  } else {
    cor_mat <- stats::cor(data, use = "pairwise.complete.obs", method = method)
  }
  return(cor_mat)
}

#' Determine the number of factors using Parallel Analysis.
#' @keywords internal
determine_n_factors <- function(data, cor_method = "pearson") {
  cor_mat <- cor_matrix_custom(data, cor_method)
  n_obs <- nrow(data)
  fa_parallel <- suppressMessages(psych::fa.parallel(
    cor_mat, fa = "fa", n.iter = 20, n.obs = n_obs,
    show.legend = FALSE, main = NULL
  ))
  return(fa_parallel$nfact)
}

#' Run a custom EFA.
#' @keywords internal
efa_custom <- function(data, n_factors = 1, cor_method = "polychoric", extract = "uls", rotate = "oblimin") {
  cor_mat <- cor_matrix_custom(data, cor_method)
  efa <- NULL
  invisible(utils::capture.output({
    efa <- psych::fa(r = cor_mat, nfactors = n_factors, rotate = rotate, fm = extract)
  }))
  alpha_val <- psych::alpha(data, check.keys = FALSE)$total$raw_alpha
  loading <- efa$loadings[]
  explained_var <- sum(efa$Vaccounted["SS loadings", 1:n_factors]) / ncol(data)
  list(efa = efa, alpha = alpha_val, explained_var = explained_var, loadings = loading)
}

#' Sort item IDs numerically.
#' @keywords internal
sort_item_ids <- function(x) {
  gtools::mixedsort(x)
}

# ' Identify problematic items (low-loading or cross-loading).
#' @keywords internal
identify_problem_items <- function(efa_res, primary_cutoff = 0.40, secondary_max = 0.30, diff_min = 0.20, low_loading_thresh = 0.30) {
  loadings <- abs(efa_res$loadings)
  items <- rownames(loadings)
  cross_3 <- c()
  cross_2 <- c()
  low_loading <- c()
  if (ncol(loadings) == 0) return(list(cross_3 = c(), cross_2 = c(), low_loading = c()))

  for (i in 1:nrow(loadings)) {
    item_loads <- loadings[i, ]
    primary <- max(item_loads, na.rm = TRUE)

    if (all(item_loads < low_loading_thresh, na.rm = TRUE)) {
      low_loading <- c(low_loading, items[i])
      next
    }

    if (ncol(loadings) > 1) {
      sorted_loads <- sort(item_loads, decreasing = TRUE)
      secondary <- sorted_loads[2]

      secondary_check <- !(primary >= primary_cutoff && secondary < secondary_max && (primary - secondary) >= diff_min)

      if (secondary_check) {
        if (sum(item_loads >= secondary_max, na.rm = TRUE) >= 3) {
          cross_3 <- c(cross_3, items[i])
        } else if (sum(item_loads >= secondary_max, na.rm = TRUE) == 2) {
          cross_2 <- c(cross_2, items[i])
        }
      }
    }
  }
  list(cross_3 = unique(cross_3), cross_2 = unique(cross_2), low_loading = unique(low_loading))
}


#' Generate combinations of items to remove.
#' @keywords internal
get_combinations <- function(items) {
  n <- length(items)
  if (n == 0) return(list(list()))
  all_combs <- unlist(lapply(1:n, function(i) utils::combn(items, i, simplify = FALSE)), recursive = FALSE)
  return(c(list(list()), all_combs)) # Add empty set for the baseline (no removal) case
}

#' Test removal strategies and build a summary table.
#' @keywords internal
test_removals <- function(data, base_items, combs, n_factors, cor_method, extract, rotate) {
  summary_list <- list()
  total <- length(combs)
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  for (i in seq_along(combs)) {
    items_to_remove <- combs[[i]]
    remaining_items <- setdiff(base_items, items_to_remove)
    if (length(remaining_items) > n_factors) {
      efa_out <- efa_custom(data[, remaining_items, drop = FALSE], n_factors, cor_method, extract, rotate)
      prob_items_test <- identify_problem_items(efa_out)

      has_cross_loading <- length(prob_items_test$cross_2) > 0 || length(prob_items_test$cross_3) > 0

      removed_label <- if (length(items_to_remove) == 0) "None" else paste(sort_item_ids(items_to_remove), collapse = "-")
      load_mat <- as.matrix(efa_out$loadings)
      threshold <- 0.30
      valid_loads <- abs(load_mat)[abs(load_mat) > threshold]
      load_min <- if (length(valid_loads) > 0) formatC(min(valid_loads), digits = 2, format = "f") else NA
      load_max <- if (length(valid_loads) > 0) formatC(max(valid_loads), digits = 2, format = "f") else NA
      loading_range <- if(is.na(load_min)) "N/A" else paste0(load_min, "-", load_max)
      summary_list[[i]] <- data.frame(
        Removed_Items = removed_label, Total_Explained_Var = efa_out$explained_var,
        Factor_Loading_Range = loading_range, Cronbachs_Alpha = efa_out$alpha,
        Cross_Loading = ifelse(has_cross_loading, "Yes", "No"), stringsAsFactors = FALSE
      )
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  summary_table <- do.call(rbind, summary_list)
  summary_table <- summary_table[!duplicated(summary_table$Removed_Items), ]

  if (nrow(summary_table) > 0) {
    sorted_indices <- order(summary_table$Cross_Loading, -summary_table$Total_Explained_Var)
    summary_table <- summary_table[sorted_indices, ]
  }

  return(summary_table)
}


# ===================================================================
# MAIN (COMPUTATION) FUNCTION
# ===================================================================

#' Evaluate Item Removal Strategies for Exploratory Factor Analysis (EFA)
#'
#' @description
#' This function automates the process of identifying low-quality items (those
#' with low factor loadings or significant cross-loadings) in an Exploratory
#' Factor Analysis (EFA). It systematically tests various combinations of
#' removing these problematic items and evaluates the impact on model fit,
#' returning a comprehensive summary of all tested strategies.
#'
#' @param data A numeric `data.frame` or `matrix` for the analysis.
#' @param cor_method The correlation method to use, e.g., `"pearson"` or `"polychoric"`.
#' @param n_factors The number of factors. If `NULL`, it's determined automatically by parallel analysis.
#' @param extract The factor extraction (estimation) method. See `psych::fa`. Default is `"uls"`.
#' @param rotate The rotation method. See `psych::fa`. Default is `"oblimin"`.
#'
#' @return An object of class `itemrest_result`. This is a list containing the
#'   following components:
#'   \item{descriptive_stats}{Basic descriptive statistics of the input data.}
#'   \item{initial_efa}{The results of the initial EFA before any items are removed.}
#'   \item{problem_items}{A list of items identified as low-loading or cross-loading.}
#'   \item{removal_summary}{A data.frame summarizing the results of all tested removal strategies.}
#'   \item{optimal_strategy}{The best-performing strategy that resulted in a clean factor structure (no cross-loadings).}
#'   \item{settings}{A list of the settings used for the analysis.}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # We will use the 'bfi' dataset from the 'psych' package.
#' # This requires the 'psych' package to be installed.
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   data(bfi, package = "psych")
#'
#'   # 1. Prepare the data: Select the personality items (first 25 columns)
#'   #    and remove rows with missing values for this example.
#'   example_data <- bfi[, 1:25]
#'   example_data <- na.omit(example_data)
#'
#'   # 2. Run the item removal analysis.
#'   #    Based on theory, the Big Five model has 5 factors.
#'   results <- itemrest(
#'     data = example_data,
#'     n_factors = 5,
#'     cor_method = "pearson" # Data is not ordinal, so pearson is appropriate
#'   )
#'
#'   # 3. Print the report for optimal strategies (default).
#'   print(results)
#'
#'   # 4. Print the report for all tested strategies.
#'   print(results, report = "all")
#' }
#' }
itemrest <- function(data,
                     cor_method = "pearson",
                     n_factors = NULL,
                     extract = "uls",
                     rotate = "oblimin") {

  # --- Parameter Validation and Setup ---
  colnames(data) <- gsub("\\.", "_", colnames(data))
  cor_method <- tolower(cor_method)
  # ... (other parameter validations would go here) ...

  n_factors_determined <- n_factors
  auto_n_factors_flag <- is.null(n_factors)
  if (auto_n_factors_flag) {
    n_factors_determined <- determine_n_factors(data, cor_method)
  }
  descriptives <- descriptive_stats(data)
  initial_efa <- efa_custom(data, n_factors_determined, cor_method, extract, rotate)

  problem_items <- identify_problem_items(initial_efa)
  all_problem_items <- sort_item_ids(unique(c(problem_items$cross_3, problem_items$cross_2, problem_items$low_loading)))

  message("--- Settings and Descriptive Statistics ---")
  if (auto_n_factors_flag) {
    message("Number of Factors (Auto): ", n_factors_determined, " (Determined by Parallel Analysis)")
  } else {
    message("Number of Factors (Manual): ", n_factors_determined)
  }
  message("Number of Items: ", descriptives$n_items)
  message("Number of Observations: ", descriptives$n_obs)
  message("Minimum Value: ", descriptives$min_value)
  message("Maximum Value: ", descriptives$max_value)
  message("\n--- Initial EFA Results (No items removed) ---")
  message("Cronbach's Alpha: ", formatC(initial_efa$alpha, digits = 3, format = "f"))
  message("Total Explained Variance: % ", formatC(initial_efa$explained_var * 100, digits = 2, format = "f"))
  message("Low-loading Items: ", ifelse(length(problem_items$low_loading) == 0, "None", paste(problem_items$low_loading, collapse = ", ")))

  cross_items_combined <- c(problem_items$cross_2, problem_items$cross_3)
  message("Cross-loading Items: ", ifelse(length(cross_items_combined) == 0, "None", paste(unique(sort_item_ids(cross_items_combined)), collapse = ", ")))

  message("\nAll Identified Low-Quality Items: ", ifelse(length(all_problem_items) == 0, "None", paste(all_problem_items, collapse = ", ")))

  # --- Test Removal Strategies ---
  removal_summary <- NULL
  optimal_strategy <- NULL

  if (length(all_problem_items) > 0) {
    all_combs <- get_combinations(all_problem_items)

    # --- Display the progress bar ---
    message("\n[Info] Testing ", length(all_combs) - 1, " different removal combinations for low-quality items...")
    removal_summary <- test_removals(data, colnames(data), all_combs, n_factors_determined, cor_method, extract, rotate)

    optimal_candidates <- removal_summary[ removal_summary[["Cross_Loading"]] == "No",
                                           , drop = FALSE ]

    if (nrow(optimal_candidates) > 0) {
      optimal_strategy <- optimal_candidates[1, , drop = FALSE]
    } else {
      optimal_strategy <- NULL
    }
  }

  # --- Create the output object ---
  output <- list(
    descriptive_stats = descriptives,
    initial_efa = initial_efa,
    problem_items = problem_items,
    all_problem_items_combined = all_problem_items,
    removal_summary = removal_summary,
    optimal_strategy = if (exists("optimal_strategy")) optimal_strategy else NULL,
    settings = list(
      n_factors = n_factors_determined, cor_method = cor_method,
      extract = extract, rotate = rotate,
      auto_n_factors = auto_n_factors_flag
    )
  )

  # Assign the custom class
  class(output) <- "itemrest_result"

  # Return the object
  return(output)
}


# ===================================================================
# PRESENTATION (PRINT) METHOD
# ===================================================================

#' Print method for `itemrest_result` class
#'
#' @param x An object of class `itemrest_result`.
#' @param report The type of report to generate: `"optimal"` (default) or `"all"`.
#' @param ... Other arguments (not used).
#'
#' @return No return value, called for side effects (prints the report to the console).
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   data(bfi, package = "psych")
#'   example_data <- na.omit(bfi[, 1:25])
#'   results <- itemrest(example_data, n_factors = 5)
#'
#'   # Print the default optimal report
#'   print(results)
#'
#'   # Print the report of all strategies
#'   print(results, report = "all")
#' }
#' }
print.itemrest_result <- function(x, report = "optimal", ...) {

  # Initial info is printed by itemrest() during computation.
  # This function only prints the final result tables and the final message.

  # Header
  cat("\n==============================\n")
  cat(" Item Removal Strategy Report\n")
  cat("==============================\n")

  if (is.null(x$removal_summary)) {
    cat("\n>> Result: No low-quality items were detected. No additional report table available.\n")
  } else {
    report_table <- NULL
    header <- ""

    # Choose which table to display based on the 'report' argument
    if (report == "all") {
      header <- "\n--- Results of All Removal Strategies ---\n"
      report_table <- x$removal_summary
    } else { # default to optimal
      header <- "\n--- Optimal Removal Strategies (No Cross-Loadings) ---\n"

      # For the optimal report, we only show strategies with No cross-loading
      report_table <- x$removal_summary[x$removal_summary$Cross_Loading == "No", ]
    }

    cat(header)
    if (is.null(report_table) || nrow(report_table) == 0) {
      cat("No strategy matching this criterion was found.\n")
    } else {

      # Make the table more readable before printing
      display_table <- report_table
      display_table$Total_Explained_Var <- paste0("% ", formatC(display_table$Total_Explained_Var * 100, format = "f", digits = 2))
      display_table$Cronbachs_Alpha <- round(display_table$Cronbachs_Alpha, 3)
      row.names(display_table) <- NULL
      print(display_table)
    }
  }

  # Final Message
  cat("\n-----------------------------------------------------\n")
  message("Final Reminder: Let algorithms be your compass, not your captain. Valid item removal also requires theoretical competence.")

  invisible(x)
}
