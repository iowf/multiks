#' Kolmogorov-Smirnov Tests for Two or More Equally-Sized Samples
#'
#' @param x a list of numeric vectors, each containing the data values for one
#'   sample
#' @param ... further arguments for methods
#'
#' @return
#' @export
#'
#' @examples
ks.rsample <- function(x, ...) {
  UseMethod("ks.rsample")
}

#' @rdname ks.rsample
#' @exportS3Method multiks::ks.rsample default
ks.rsample.default <- function(
  x,
  exact = NULL,
  simulate.p.value = FALSE,
  B = 2e3
) {

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

  k <- maxdist(x, r)
  D <- k / n

  method <-
    if ((r * k)^r < 1e4) "Exact"
    else if (simulate.p.value) "Monte-Carlo"
    else "Approximate"

  if (method == "Approximate" && n < 10) warning("Approximate p-value may be innaccurate for small sample sizes")

  P <- switch(method,
    "approximate" = p_ks_approx(r, n, k),
    "exact" = p_ks_exact(r, n, k),
    "simulated" = p_ks_sim(r, n, k, B)
  )

  structure(
    list(
      statistic = D,
      p.value = P,
      method = paste(method, "r-sample Kolmogorov-Smirnov test")
    ),
    class = c("ks.rsample", "ks.test", "htest")
  )
}

maxdist <- function(x) {
  stopifnot(is.list(x))
  stopifnot(all(vapply(x, is.numeric, logical(1L))))
  r <- length(x)
  N <- vapply(x, length, integer(1L))
  n <- unique(N)
  stopifnot(length(n) == 1L)

  samples <- rep(1L:r, each = n) |> unlist()
  ascending <- samples[order(unlist(x))]
  lattice <- outer(ascending, 1L:r, `==`) |> apply(2L, cumsum)
  circdiffs <- lapply(1L:r, \(i) lattice[, i] - lattice[, i %% r + 1L]) |>
    do.call(cbind, args = _)
  max(circdiffs)
}
