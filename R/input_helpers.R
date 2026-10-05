# Internal helpers shared by search, validation, and resampling.
check_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " must be a single ", if (integer) "integer" else "number",
         " between ", lower, " and ", upper, ".", call. = FALSE)
  }
  invisible(x)
}

resolve_itemrest_seed <- function(seed) {
  if (is.null(seed)) {
    label <- format(Sys.Date(), "%d%m%Y")
    return(list(value = as.integer(label), label = label, source = "calendar_date"))
  }
  check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  list(value = as.integer(seed), label = as.character(seed), source = "user")
}

apply_item_keys <- function(data, keys = NULL) {
  full <- stats::setNames(rep(1, ncol(data)), names(data))
  if (!is.null(keys)) {
    if (!is.numeric(keys) || is.null(names(keys)) || anyNA(keys) || anyDuplicated(names(keys)) ||
        any(!keys %in% c(-1, 1)) || any(!names(keys) %in% names(data)))
      stop("keys must be a named numeric vector of 1 or -1 for known items.", call. = FALSE)
    full[names(keys)] <- keys
  }
  for (item in names(data)) data[[item]] <- data[[item]] * full[[item]]
  list(data = data, keys = full)
}

with_itemrest_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

capture_analysis <- function(code) {
  notes <- character()
  messages <- character()
  error <- NULL
  error_condition <- NULL
  value <- withCallingHandlers(
    tryCatch(force(code), error = function(e) {
      error <<- conditionMessage(e)
      error_condition <<- e
      NULL
    }),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, warnings = unique(notes), messages = unique(messages),
       error = error, error_condition = error_condition)
}

validate_efa_methods <- function(extract, rotate) {
  methods <- c("minres", "uls", "ols", "wls", "gls", "pa", "ml", "minchi", "minrank", "alpha", "old.min")
  rotations <- c("none", "varimax", "quartimax", "bentlerT", "equamax", "varimin",
    "geominT", "bifactor", "promax", "Promax", "oblimin", "simplimax", "bentlerQ", "geominQ",
    "biquartimin", "cluster")
  resolve <- function(value, choices, name) {
    if (!is.character(value) || length(value) != 1L || is.na(value))
      stop(name, " must be one supported method name.", call. = FALSE)
    if (value %in% choices) return(value)
    index <- match(tolower(value), tolower(choices))
    if (is.na(index)) stop("Unsupported ", name, " = '", value,
      "'. Choose from: ", paste(choices, collapse = ", "), ".", call. = FALSE)
    choices[index]
  }
  extract <- resolve(extract, methods, "extract")
  rotate <- resolve(rotate, rotations, "rotate")
  has_gpa <- rotate != "none" && requireNamespace("GPArotation", quietly = TRUE)
  if (!rotate %in% c("none", "varimax", "promax", "Promax", "cluster") && !has_gpa)
    stop("rotate = '", rotate, "' requires the GPArotation package.", call. = FALSE)
  list(extract = extract, rotate = rotate)
}

# Informational messages are recorded without changing candidate eligibility.
# Only these explicit estimator problem reports require review.
estimation_message_flags <- function(messages) {
  patterns <- c(
    nonconvergence = "failed.*converg|not.*converg|converg.*(fail|not|limit)|iteration limit|maximum.*iter",
    rotation_unavailable = "rotation.*(not found|not available|unsupported)|specified rotation not found",
    reported_variance_problem = "heywood|negative.*(uniqueness|variance)|non.?positive.*uniqueness",
    reported_matrix_repair = "(matrix|correlation).*(smooth|repair)|smoothing.*(done|applied)"
  )
  names(patterns)[vapply(patterns, function(pattern)
    any(grepl(pattern, messages, ignore.case = TRUE)), logical(1))]
}

prepare_itemrest_data <- function(data, missing = "listwise") {
  if (!(is.data.frame(data) || is.matrix(data)))
    stop("data must be a numeric data.frame or matrix.", call. = FALSE)
  if (is.matrix(data) && !is.numeric(data))
    stop("All items must be numeric.", call. = FALSE)
  data <- as.data.frame(data, check.names = FALSE)
  if (ncol(data) < 2L) stop("data must contain at least two items.", call. = FALSE)
  nonnumeric <- names(data)[!vapply(data, is.numeric, logical(1))]
  if (length(nonnumeric)) stop("All items must be numeric. Nonnumeric items: ",
    paste(nonnumeric, collapse = ", "), ".", call. = FALSE)
  item_names <- names(data)
  if (anyNA(item_names) || any(!nzchar(item_names)) || anyDuplicated(item_names))
    stop("Item names must be nonempty and unique.", call. = FALSE)
  if (any(vapply(data, function(x) any(is.infinite(x)), logical(1))))
    stop("Infinite item values are not supported.", call. = FALSE)
  missing <- match.arg(missing, c("listwise", "pairwise", "fail"))
  n_input <- nrow(data)
  n_missing <- sum(!stats::complete.cases(data))
  if (missing == "fail" && n_missing > 0L)
    stop("Missing values found with missing = 'fail'.", call. = FALSE)
  rows <- if (missing == "listwise") which(stats::complete.cases(data)) else seq_len(n_input)
  data <- data[rows, , drop = FALSE]
  if (nrow(data) < 3L) stop("At least three observations are required after missing-data handling.", call. = FALSE)
  list(data = data, summary = list(
    n_input = n_input, n_used = nrow(data), n_excluded = n_input - nrow(data),
    n_rows_with_missing = n_missing, missing = missing, retained_rows = rows
  ))
}

pairwise_counts <- function(data) {
  observed <- !is.na(as.matrix(data))
  crossprod(observed * 1L)
}

is_ordinal_item <- function(x, ordinal_categories = 7L) {
  x <- x[!is.na(x)]
  length(unique(x)) <= ordinal_categories && length(x) > 0L && all(x %% 1 == 0)
}

item_label <- function(items) {
  if (length(items) == 0L) "None" else paste(gtools::mixedsort(items), collapse = "-")
}

empty_problem_items <- function() {
  list(cross_3 = character(), cross_2 = character(), low_loading = character())
}

all_flagged_items <- function(problems) {
  unique(c(problems$cross_3, problems$cross_2, problems$low_loading))
}

analysis_options <- function(missing, min_items_per_factor, max_factor_correlation,
                             pd_action, reliability) {
  list(missing = missing, min_items_per_factor = min_items_per_factor,
       max_factor_correlation = max_factor_correlation, pd_action = pd_action,
       reliability = reliability)
}

package_provenance <- function() {
  packages <- c("ItemRest", "psych", "qgraph", "gtools")
  versions <- vapply(packages, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  list(R = R.version.string, packages = versions, RNG_kind = RNGkind(),
       timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

data_fingerprint <- function(data, ignore_row_order = FALSE) {
  values <- unname(as.matrix(data))
  storage.mode(values) <- "double"
  values[is.na(values)] <- NA_real_
  values[!is.na(values) & values == 0] <- 0
  if (ignore_row_order) {
    index <- do.call(order, c(unname(as.list(as.data.frame(values))),
      list(na.last = TRUE, method = "radix")))
    values <- values[index, , drop = FALSE]
  }
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(items = names(data), values = values), path, version = 2)
  unname(tools::md5sum(path))
}
