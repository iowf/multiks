
p_ks_exact <- function(r, n, k) {
  # term out front
  comb <- (factorial(n)^r) / ((r * k)^r * factorial(n * r))
  # each row is one combination of v1, ..., vr
  vv <- expand.grid(
    rep(list(1:(r*k) - 1), r)
  ) |> as.matrix()
  # remove unneeded indices
  upper <- upper.tri(diag(rep(1, r)))
  ii <- outer(1:r, rep(1, r))[upper]
  jj <- outer(rep(1, r), 1:r)[upper]
  diffs <- mapply(
    \(i, j) vv[, i] - vv[, j],
    ii, jj,
    SIMPLIFY = TRUE
  )
  becomes_one <- apply(diffs %% r == 0, 1, any)
  vv <- vv[!becomes_one, ]
  # omega
  w <- exp((2 * pi * 1i) / (r * k))
  # first term within sum
  sum2 <- rowSums(w ^ vv)^(n * r)
  # second term within sum
  exp1 <- w^(-n * rowSums(vv))
  # third term within sum
  prod1 <- apply(vv, 1, \(v) prod(
    1 - w^(k*outer(v, v, `-`)[upper])
  ))
  # full expression
  P <- comb * sum(sum2 * exp1 * prod1)
  stopifnot(isTRUE(all.equal(Im(P), 0)))
  Re(P)
}

p_ks_approx <- function(r, n, k) {
  diffs <- {mat <- outer(1:r, 1:r, `-`); mat[lower.tri(mat)]}
  factorial(n)^r / factorial(r * n) *
    2 ^ (r * (r - 1)) / (r * k)^(r - 1) *
    (sin(pi / k) / sin(pi / (r * k)))^(r * n) *
    prod(sin(pi * diffs / r)^2)
}
