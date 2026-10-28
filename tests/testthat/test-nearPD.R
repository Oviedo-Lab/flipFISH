# File written by Claude Sonnet 5

# Regression coverage for dg_nearPD_native() -- the native C++ port of R's Matrix::nearPD()
# (Higham 2002's alternating-projections algorithm with Dykstra's correction) that
# dg_nearest_corr_matrix() now uses in place of an R callback (see src/main.cpp). This checks that the
# port agrees with Matrix::nearPD(corr = TRUE, keepDiag = TRUE, do2eigen = TRUE) -- the exact call this
# package has always made -- across matrices shaped like what it's actually used on: random
# indefinite symmetric matrices, matrices built the same pairwise-independent way
# dg_raw_gaussian_corr() builds them (unit diagonal, each off-diagonal entry chosen without any joint
# PSD guarantee), and already-valid correlation matrices (which should pass through essentially
# unchanged).

skip_if_not_installed("Matrix")

tol <- 1e-8 # observed agreement during development was ~1e-14 (float roundoff); 1e-8 leaves margin

check_against_nearPD <- function(x) {
  ref <- as.matrix(Matrix::nearPD(x, corr = TRUE, keepDiag = TRUE, do2eigen = TRUE)$mat)
  got <- flipFISH:::dg_nearPD_native_test(x)
  expect_equal(got, ref, tolerance = tol)
}

random_indefinite <- function(n, seed) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n, n)
  S <- (A + t(A)) / 2
  diag(S) <- 1
  S
}

pairwise_style <- function(n, seed) {
  set.seed(seed)
  S <- diag(n)
  for (i in 2:n) {
    for (j in seq_len(i - 1)) {
      rho <- runif(1, -0.95, 0.95)
      S[i, j] <- rho
      S[j, i] <- rho
    }
  }
  S
}

test_that("dg_nearPD_native matches Matrix::nearPD on random indefinite matrices", {
  for (n in c(2, 3, 5, 10, 20, 28)) {
    for (seed in 1:3) {
      check_against_nearPD(random_indefinite(n, seed))
    }
  }
})

test_that("dg_nearPD_native matches Matrix::nearPD on dg_raw_gaussian_corr()-shaped matrices", {
  # This is the actual input shape in production: unit diagonal, each off-diagonal entry solved
  # independently per bit-pair (dg_raw_gaussian_corr()), so the assembled matrix has no joint PSD
  # guarantee -- exactly the case dg_nearest_corr_matrix()'s projection step exists to handle.
  for (n in c(4, 8, 16, 28)) {
    for (seed in 1:3) {
      check_against_nearPD(pairwise_style(n, seed))
    }
  }
})

test_that("dg_nearPD_native leaves an already-valid correlation matrix essentially unchanged", {
  set.seed(1)
  for (n in c(3, 8, 28)) {
    A <- matrix(rnorm(n * n), n, n)
    S <- cov2cor(A %*% t(A) + diag(n) * 0.1)
    check_against_nearPD(S)
  }
})
