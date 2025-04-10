test_that("calculated values match authors' values", {
  p_table <- read.delim("data/exactapprox.tsv") # Table 1 from Böhm & Hornik
  expected <- p_table$P.nd...k.
  funs <- list("exact" = p_rks_exact, "approx" = p_rks_approx)
  actual <- split(p_table, seq_len(nrow(p_table))) |> vapply(\(row) {
    with(row, funs[[method]](r, n, k))
  }, numeric(1), USE.NAMES = FALSE)
  expect_equal(
    expected,
    round(actual, 4)
  )
})

test_that("gives same results as ks.test() for 2 samples", {
  set.seed(5687)
  x <- rnorm(15)
  y <- runif(15)
  expect_equal(
    ks_rsample(list(x, y))[c("statistic", "p.value")],
    ks.test(x, y)[c("statistic", "p.value")]
  )
})

test_that("input validation works", {
  expect_error(regexp = "size",
    ks_rsample(list(runif(10), runif(15)))
  )
  expect_error(regexp = "ties",
    ks_rsample(list(c(1, 2, 4, 7, 12), c(3, 5, 7, 8, 10)))
  )
})

test_that("methods are equivalent", {
  test_methods <- function(x) {
    x_list <- x
    x_matrix <- do.call(cbind, args = x_list)
    x_data.frame <- x_list |>
      as.data.frame() |>
      stats::setNames(1:4) |>
      utils::stack()

    ks_matrix <- ks_rsample(x_matrix)
    ks_matrix$data.name <- NULL
    ks_list <- ks_rsample(x_list)
    ks_list$data.name <- NULL
    ks_data.frame <- ks_rsample(values ~ ind, data = x_data.frame)
    ks_data.frame$data.name <- NULL

    expect_equal(ks_matrix, ks_list)
    expect_equal(ks_list, ks_data.frame)
    expect_equal(ks_data.frame, ks_matrix)
  }

  set.seed(8675309)
  # for exact p-values
  x1 <- lapply(1:4, \(i) runif(5))
  test_methods(x1)
  # for inexact p-values
  x2 <- lapply(1:4, \(i) runif(25))
  test_methods(x2)
})

# create a weird distribution to capture some edge cases
weird_dist <- function(k = 5) {
  force(k)
  means <- 2 * rnorm(k)
  sds <- exp(rnorm(k))
  weights <- runif(k)
  weights <- weights / sum(weights)
  function(n) {
    which_bells <- sample(1:k, n, replace = TRUE, prob = weights)
    rnorm(n, means[which_bells], sds[which_bells])
  }
}
test_that("simulated p-values are reasonably accurate and precise", {
  set.seed(03755)
  iter <- 100
  rr <- pmax(floor(10^runif(iter, 0, 1)), 2)
  nn <- ceiling(10^runif(iter, 1, 2))
  p_values <- mapply(rr, nn, FUN = \(r, n) {
    x <- lapply(1:r, \(i2) (weird_dist())(n))
    tryCatch(
      {
        calculated <- ks_rsample(x, exact = NULL, simulate.p.value = FALSE)
        simulated  <- ks_rsample(x, exact = FALSE, simulate.p.value = TRUE)
        c(calculated = calculated$p.value, simulated = simulated$p.value)
      },
      warning = \(w) {
        stopifnot(grepl("inaccurate", w$message))
        return(NULL)
      }
    )
  }, SIMPLIFY = FALSE) |> do.call(rbind, args = _) |> as.data.frame()
  diffs <- with(p_values, (simulated - calculated) / calculated)
  # expect that simulated p-values are not biased
  expect_gte(t.test(diffs)$p.value, 0.05)
  # expect that simulated p-values are within +/- 10% of the calculated p-value
  #   90% of the time (a low bar, but the default number of iterations is
  #   fairly low...)
  expect_gte(mean(abs(diffs) < 0.1), 0.9)
})
