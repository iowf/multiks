# Interfaces ===================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Kolmogorov-Smirnov Tests for Two or More Equally-Sized Samples
#'
#' An implementation of Walter Böhm and Kurt Hornik's (2010) two-or-more-sample
#' Kolmogorov-Smirnov test
#'
#' @param x a list of numeric vectors, each containing the data values for one
#'   sample; or a matrix, with each column containing the data values for one
#'   sample. Each sample must be the same size, and the data must not include
#'   missing vaues or ties.
#' @param formula a formula, where the left-hand side is a numeric variable in
#'   `data` and the right-hand side contains a single term whose unique values
#'   define the groups for the test. The LHS must not contain any missing
#'   [model.frame()].
#' @param data an optional data frame (or similar; see [model.frame()])
#'   containing variables used in `formula`.
#' @param exact logical; if `TRUE`, or if `NULL` and `(r * D * n) ^ r < 1e4`,
#'   then an exact p-value is computed; otherwise, an inexact p-value is
#'   computed according to the setting of `simulate.p.value`. (`r` is the number
#'   of samples, `D` is the test statistic, and `n` is the size of each sample)
#' @param simulate.p.value logical; ignored in an exact p-value is computed. If
#'   `TRUE`, the p-value is estimated using Monte Carlo simulation, and if
#'   `FALSE`, using the approximation from Böhm and Hornik (2011). `FALSE` is
#'   recommended except for small sample sizes, as the approximation is fairly
#'   good at least down to sample sizes of about 10.
#' @param B integer; number of trials for the Monte Carlo p-value estimation
#' @param ... further arguments for methods
#'
#' @return An object of class "ks_rsample".
#'
#' @details For the purpose of calculating the test statistic, `model.formula()`
#'   orders the groups by calling [order()] on the values of the grouping
#'   variables. If the right-hand side of the formula is an interaction term,
#'   each unique level of the interaction is between those variables is taken to
#'   be one group.
#'
#' @references Böhm, W., & Hornik, K. (2010). A Kolmogorov-Smirnov Test for r
#'   Samples. \emph{Research Report Series / Department of Statistics and
#'   Mathematics} 105.
#'   \href{https://doi.org/10.57938/2d4c243c-f18c-450e-8d1d-d0a69ff1fe52}{doi:10.57938/2d4c243c-f18c-450e-8d1d-d0a69ff1fe52}
#'
#' @export
#'
#' @examples
#' x <- list(rnorm(100, sd = 1), rnorm(100, sd = 2), rnorm(100, sd = 3))
#' ks_rsample(x)
#'
#' # using a formula and data frame:
#' ks_rsample(Sepal.Length ~ Species, data = iris)
ks_rsample <- function(x, ...) {
  UseMethod("ks_rsample")
}

#' @rdname ks_rsample
#' @exportS3Method multiks::ks_rsample default
ks_rsample.default <- function(
  x,
  exact = NULL, simulate.p.value = NULL, B = 1e4,
  ...
) {

  data_name <- switch(match("data_name", ...names()), ...)
  if (is.null(data_name)) data_name <- deparse1(substitute(x))

  # input validation
  if (!is.list(x)) stop(sprintf("'x' must be a list, not %s", typeof(x)), call. = FALSE)
  if (length(x) < 2) stop(sprintf("'x' has %d elements, but requires at least 2", length(x)), call. = FALSE)
  if (!all(vapply(x, is.numeric, logical(1)))) stop("All elements of 'x' must be numeric", call. = FALSE)
  if (any(is.na(unlist(x)))) stop("'x' contains NA values", call. = FALSE)
  if (any(duplicated(unlist(x)))) stop("'x' contains ties", call. = FALSE)

  r <- length(x)
  n <- vapply(x, length, integer(1)) |> unique()

  if (length(n) > 1) stop("All samples must be of the same size", call. = FALSE)
  if (n < 1) stop("Not enough 'x' values", call. = FALSE)

  k <- max_dist(x)
  D <- k / n

  # check limitations on methods for computing P
  exact_iter_too_big <- (r * k)^r >= 1e6 # will take very long to compute
  approx_k_too_big <- k > sqrt(n) # inaccurate for large k (anecdotal)
  approx_n_too_small <- n < 10 # not tested for small n

  method <-
    # if exact = TRUE, always use 'exact'
    if (isTRUE(exact)) "Exact"
    # if exact is NULL, then use 'exact' if practical
    else if (!isFALSE(exact) && !exact_iter_too_big) "Exact"
    # if simulate.p.value = FALSE, always use 'approx'
    else if (isFALSE(simulate.p.value)) "Approximate"
    # if simulate.p.value is NULL, use 'approx' if reasonably accurate
    else if (!isTRUE(simulate.p.value) &&
      !approx_k_too_big && !approx_n_too_small) "Approximate"
    # if specified, or if all else fails, use 'sim'
    else "Monte-Carlo"

  if (method == "Approximate") {
    if (approx_k_too_big) warning("Approximate p-value may be inaccurate when ~ D > 1/sqrt(n)")
    if (approx_n_too_small) warning("Approximate p-value may be inaccurate for small sample sizes")
  }

  P <- 1 - switch(method,
    "Exact" = p_rks_exact(r, n, k),
    "Monte-Carlo" = p_rks_sim(r, n, k, B),
    "Approximate" = p_rks_approx(r, n, k),
  )

  if (method == "Monte-Carlo" && P * (B + 1) < 100) warning("Fewer than 100 iterations had greater than the observed D; simulated P-value is likely inaccurate. Try increasing the number of iterations.")

  structure(
    list(
      statistic = c(D = D),
      p.value = P,
      method = paste(method, "r-sample Kolmogorov-Smirnov test"),
      data.name = data_name,
      exact = (method == "Exact"),
      iterations = if (method == "Monte-Carlo") B else NA
    ),
    class = c("ks_rsample", "ks.test", "htest")
  )
}

#' @rdname ks_rsample
#' @exportS3Method multiks::ks_rsample matrix
ks_rsample.matrix <- function(
  x,
  exact = NULL, simulate.p.value = FALSE, B = 1e4,
  ...
) {
  data_name <- deparse1(substitute(x))
  if (!is.numeric(x)) stop(sprintf("'x' must be numeric, not %s", typeof(x)), call. = FALSE)
  if (is.matrix(x) & ncol(x) < 2) stop(sprintf("'x' has %d columns, but requires at least 2", ncol(x)), call. = FALSE)
  x <- x |> apply(2, identity, simplify = FALSE)
  NextMethod(x = x, data_name = data_name, ...)
}

#' @rdname ks_rsample
#' @exportS3Method multiks::ks_rsample formula
ks_rsample.formula <- function(
  formula, data,
  exact = NULL, simulate.p.value = FALSE, B = 1e4,
  ...
){
  # browser()
  mframe <- stats::model.frame(formula, data = data)
  mterms <- attributes(stats::terms(mframe))
  if (ncol(mterms$factors) != 1) stop(sprintf("'formula' must contain 1 term on RHS, not %d", ncol(mterms$facotrs)), call. = FALSE)
  if (!is.numeric(mframe[[mterms$response]])) stop(sprintf("'%s' must be numeric, not %s", as.character(mterms$variables[[mterms$response + 1]]), typeof(mframe[[mterms$response]])), call. = FALSE)

  mframe <- mframe[order(mframe[, mterms$factors[, 1] > 0]), ]
  x <- split(mframe[[mterms$response]], mframe[, mterms$factors[, 1] > 0])

  data_name <- deparse1(formula)

  NextMethod(x = x, data_name = data_name, ...)
}

# Auxiliary functions ==========================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate k, n times the maximum circular difference between ECDFs of the
#   samples
#   x   A list of equally-sized numeric vectors
max_dist <- function(x) {
  stopifnot(is.list(x))
  stopifnot(all(vapply(x, is.numeric, logical(1L))))
  r <- length(x)
  n <- vapply(x, length, integer(1L)) |> unique()
  stopifnot(length(n) == 1L)

  samp_nums <- rep(1L:r, each = n) |> unlist()
  ascending <- samp_nums[order(unlist(x))]
  max_dist_tight(ascending, r)
}
# The workhorse of max_dist(); no overhead of validating and sorting x
#   ascending   An integer vector; the sample indices (1, ..., r) of the
#               observations in x sorted by the values of the observations in x
#   r           The number of samples
max_dist_tight <- function(ascending, r) {
  lattice <- outer(ascending, 1L:r, `==`) |> apply(2L, cumsum)
  circdiffs <- lapply(1L:r, \(i) lattice[, i] - lattice[, i %% r + 1L]) |>
    do.call(cbind, args = _)
  max(circdiffs)
}

# Probabilities ================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Each function calculates the probability that n times the maximum circular
#   difference between sample ECDFs is less than k, for r samples of size n

# exact probability (Eq. 17 from Böhm & Hornik)
p_rks_exact <- function(r, n, k) {
  # each row is one combination of the summation indices v1, ..., vr
  vv <- expand.grid(
    rep(list(1:(r*k) - 1), r)
  ) |> as.matrix()
  # remove unneeded indices
  upper <- upper.tri(diag(rep(1, r)))
  diffs <- mapply(
    \(i, j) vv[, i] - vv[, j],
    outer(1:r, rep(1, r))[upper],
    outer(rep(1, r), 1:r)[upper],
    SIMPLIFY = TRUE
  )
  becomes_one <- apply(diffs %% r == 0, 1, any)
  vv <- vv[!becomes_one, ]
  w <- exp((2 * pi * 1i) / (r * k)) # omega
  P <- exp(r * lgamma(n + 1) - r * log(r * k) - lgamma(n * r + 1)) * sum(
    rowSums(w ^ vv)^(n * r) *
    w^(-n * rowSums(vv)) *
    apply(vv, 1, \(v) prod(
      1 - w^(k*outer(v, v, `-`)[upper])
    ))
  )
  if (!isTRUE(all.equal(Im(P), 0))) browser()
  stopifnot(isTRUE(all.equal(Im(P), 0)))
  Re(P)
}

# simulated probability (Monte-Carlo)
#   B   number of iterations
p_rks_sim <- function(r, n, k, B) {
  samp_nums <- rep(1L:r, each = n) |> unlist()
  N <- n * r
  nge <- vapply(1:B, \(i)
    max_dist_tight(sample(samp_nums, N), r) < k,
  logical(1)) |> sum()
  (nge + 1) / (B + 1)
}

# approximation (Eq. 16 from Böhm & Hornik)
p_rks_approx <- function(r, n, k) {
  diffs <- {mat <- outer(1:r, 1:r, `-`); mat[lower.tri(mat)]}
  exp(
    lgamma(n + 1)*r - lgamma(r * n + 1) +
    (r * (r - 1)) * log(2) - (r - 1) * log(r * k) +
    (r * n) * log(sin(pi / k) / sin(pi / (r * k))) +
    sum(2 * log(sin(pi * diffs / r)))
  )
}
