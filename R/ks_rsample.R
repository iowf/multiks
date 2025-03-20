#' Kolmogorov-Smirnov Tests for Two or More Equally-Sized Samples
#'
#' @param x a list of numeric vectors, each containing the data values for one
#'   sample; or a matrix, with each column containing the data values for one sample.
#'   Each sample must be the same size and must not contain NAs.
#' @param exact logical; if `TRUE`, or if `NULL` and `(r * D * n) ^ r < 1e4`,
#'   then an exact p-value is computed; otherwise, an inexact p-value is computed
#'   according to the setting of `simulate.p.value`. (`r` is the number of samples,
#'   `D` is the test statistic, and `n` is the size of each sample)
#' @param simulate.p.value logical; ignored in an exact p-value is computed. If
#'   `TRUE`, the p-value is estimated using Monte Carlo simulation, and if `FALSE`,
#'   using the approximation from Böhm and Hornik (2011). `FALSE` is recommended
#'   except for small sample sizes, as the approximation is fairly good
#'   at least down to sample
#'   sizes of about 10.
#' @param B integer; number of trials for the Monte Carlo p-value estimation
#' @param ... further arguments for methods
#'
#' @return An object of class "ks.rsample".
#' @export
#'
#' @examples
#' # using a list as input:
#' x <- list(runif(100, 0), runif(100, 0.5), rnorm(100, 1))
#' ks.rsample(x)
ks.rsample <- function(x, ...) {
  UseMethod("ks.rsample")
}

#' @rdname ks.rsample
#' @exportS3Method multiks::ks.rsample default
ks.rsample.default <- function(
  x,
  exact = NULL,
  simulate.p.value = FALSE,
  B = 2e3,
  ...
) {
  data_name <- deparse1(substitute(x))
  browser()
  # input validation
  if (is.list(x)) {
    if (length(x) < 2) stop(sprintf("'x' has %d elements, but requires at least 2", length(x)))
    if (!all(vapply(x, is.numeric, logical(1)))) stop("All elements of 'x' must be numeric")
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) stop(sprintf("'x' must be numeric, not %d", typeof(x)))
    if (is.matrix(x) & ncol(x) < 2) stop(sprintf("'x' has %d columns, but requires at least 2", ncol(x)))
    x <- x |> apply(2, identity, simplify = FALSE)
  } else stop(sprintf("'x' must be a list or matrix, not a %s", typeof(x)))
  if (any(is.na(unlist(x)))) stop("'x' contains NA values")

  r <- length(x)
  n <- vapply(x, length, integer(1)) |> unique()

  if (length(n) > 1) stop("All samples must be of the same size")
  if (n < 1) stop("Not enough 'x' values")

  k <- max_dist(x)
  D <- k / n

  method <-
    if (isTRUE(exact)) "Exact"
    else if (!isFALSE(exact) && (r * k)^r < 1e4) "Exact"
    else if (simulate.p.value) "Monte-Carlo"
    else "Approximate"
  if (method == "Approximate" && n < 10) warning("Approximate p-value may be innaccurate for small sample sizes")

  P <- switch(method,
    "Exact" = p_rks_exact(r, n, k),
    "Monte-Carlo" = p_rks_sim(r, n, k, B),
    "Approximate" = p_rks_approx(r, n, k),
  )

  structure(
    list(
      statistic = c(D = D),
      p.value = P,
      method = paste(method, "r-sample Kolmogorov-Smirnov test"),
      data.name = data_name,
      exact = (method == "Exact")
    ),
    class = c("ks.rsample", "ks.test", "htest")
  )
}

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
max_dist_tight <- function(ascending, r) {
  lattice <- outer(ascending, 1L:r, `==`) |> apply(2L, cumsum)
  circdiffs <- lapply(1L:r, \(i) lattice[, i] - lattice[, i %% r + 1L]) |>
    do.call(cbind, args = _)
  max(circdiffs)
}

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
  # omega
  w <- exp((2 * pi * 1i) / (r * k))
  # 1 - P(nD < k)
  Pc <- (factorial(n)^r) / ((r * k)^r * factorial(n * r)) * sum(
    rowSums(w ^ vv)^(n * r) *
    w^(-n * rowSums(vv)) *
    apply(vv, 1, \(v) prod(
      1 - w^(k*outer(v, v, `-`)[upper])
    ))
  )
  stopifnot(isTRUE(all.equal(Im(Pc), 0)))
  1 - Re(Pc)
}

p_rks_sim <- function(r, n, k, B) {
  samp_nums <- rep(1L:r, each = n) |> unlist()
  N <- n * r
  nge <- vapply(1:B, \(i)
    max_dist_tight(sample(samp_nums, N), r) >= k,
  logical(1)) |> sum()
  (nge + 1) / (B + 1)
}

p_rks_approx <- function(r, n, k) {
  diffs <- {mat <- outer(1:r, 1:r, `-`); mat[lower.tri(mat)]}
  exp(
    lgamma(n + 1)*r - lgamma(r * n + 1) +
    (r * (r - 1)) * log(2) - (r - 1) * log(r * k) +
    (r * n) * log(sin(pi / k) / sin(pi / (r * k))) +
    sum(2 * log(sin(pi * diffs / r)))
  )
}
