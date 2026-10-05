fixture_data <- function(n = 100L, p = 9L) {
  matrix <- outer(seq_len(n), seq_len(p), function(i, j) sin(i * j / 7))
  as.data.frame(matrix, optional = TRUE) |>
    stats::setNames(paste0("I", seq_len(p)))
}

make_test_efa <- function(items, low = character(), cross = character(), phi = NULL) {
  loadings <- matrix(0.05, nrow = length(items), ncol = 2L,
                     dimnames = list(items, c("F1", "F2")))
  primary <- rep(1:2, length.out = length(items))
  loadings[cbind(seq_along(items), primary)] <- 0.70
  loadings[items %in% low, ] <- 0.10
  if (length(cross)) loadings[items %in% cross, ] <- rep(c(0.50, 0.48), each = sum(items %in% cross))
  list(efa = list(converged = TRUE), alpha = 0.8, explained_var = 0.5,
       loadings = loadings, phi = if (is.null(phi)) diag(2) else phi,
       correlation = list(matrix = diag(length(items)), positive_definite = TRUE,
                          adjusted = FALSE, min_eigenvalue = 1), fit_warnings = character())
}

real_factor_data <- function(n = 300L, p = 5L, seed = 22L) {
  with_itemrest_seed(seed, {
    f <- rnorm(n)
    as.data.frame(replicate(p, f + rnorm(n, sd = 0.6)))
  })
}
