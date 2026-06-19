
# File written by Claude Sonnet 4.6

# Tests for the TR() probability model.
# tr_sum_check() sums TR(bc, flips, fr) over all 2^N flip patterns for a given barcode;
# the result must be 1.0 (within floating-point tolerance) for the model to be a valid
# probability distribution.

library(flipFISH)

tol <- 1e-12   # tolerance for sum-to-1 checks

# ---- helpers ----------------------------------------------------------------

# Build zero correlation vectors for N bits
no_corr <- function(N) rep(0, N * (N - 1) / 2)

# ---- 2-bit tests ------------------------------------------------------------

test_that("TR sums to 1 for 2 bits, zero correlations, bc = 00", {
  N <- 2
  expect_equal(
    tr_sum_check(bc = 0L,
                 rate10 = c(0.05, 0.08),
                 rate01 = c(0.03, 0.06),
                 corr1  = no_corr(N),
                 corr0  = no_corr(N)),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 2 bits, zero correlations, bc = 11", {
  N <- 2
  expect_equal(
    tr_sum_check(bc = 3L,
                 rate10 = c(0.05, 0.08),
                 rate01 = c(0.03, 0.06),
                 corr1  = no_corr(N),
                 corr0  = no_corr(N)),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 2 bits, non-zero correlations, bc = 00", {
  N <- 2
  expect_equal(
    tr_sum_check(bc = 0L,
                 rate10 = c(0.05, 0.08),
                 rate01 = c(0.03, 0.06),
                 corr1  = 0.3,
                 corr0  = 0.2),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 2 bits, non-zero correlations, bc = 11", {
  N <- 2
  expect_equal(
    tr_sum_check(bc = 3L,
                 rate10 = c(0.05, 0.08),
                 rate01 = c(0.03, 0.06),
                 corr1  = 0.3,
                 corr0  = 0.2),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 2 bits, non-zero correlations, bc = 01", {
  N <- 2
  expect_equal(
    tr_sum_check(bc = 1L,
                 rate10 = c(0.05, 0.08),
                 rate01 = c(0.03, 0.06),
                 corr1  = 0.3,
                 corr0  = 0.2),
    1.0, tolerance = tol
  )
})

# ---- 4-bit tests ------------------------------------------------------------

test_that("TR sums to 1 for 4 bits, zero correlations, bc = 0101", {
  N <- 4
  expect_equal(
    tr_sum_check(bc = 5L,   # 0101 in binary
                 rate10 = c(0.05, 0.08, 0.04, 0.07),
                 rate01 = c(0.03, 0.06, 0.02, 0.05),
                 corr1  = no_corr(N),
                 corr0  = no_corr(N)),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 4 bits, non-zero correlations, bc = 1111", {
  N <- 4
  # 6 lower-triangle entries for 4 bits
  expect_equal(
    tr_sum_check(bc = 15L,
                 rate10 = c(0.05, 0.08, 0.04, 0.07),
                 rate01 = c(0.03, 0.06, 0.02, 0.05),
                 corr1  = c(0.20, 0.15, 0.10, 0.25, 0.05, 0.18),
                 corr0  = c(0.10, 0.08, 0.12, 0.09, 0.15, 0.07)),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 4 bits, non-zero correlations, bc = 0000", {
  N <- 4
  expect_equal(
    tr_sum_check(bc = 0L,
                 rate10 = c(0.05, 0.08, 0.04, 0.07),
                 rate01 = c(0.03, 0.06, 0.02, 0.05),
                 corr1  = c(0.20, 0.15, 0.10, 0.25, 0.05, 0.18),
                 corr0  = c(0.10, 0.08, 0.12, 0.09, 0.15, 0.07)),
    1.0, tolerance = tol
  )
})

test_that("TR sums to 1 for 4 bits, non-zero correlations, bc = 1010", {
  N <- 4
  expect_equal(
    tr_sum_check(bc = 10L,
                 rate10 = c(0.05, 0.08, 0.04, 0.07),
                 rate01 = c(0.03, 0.06, 0.02, 0.05),
                 corr1  = c(0.20, 0.15, 0.10, 0.25, 0.05, 0.18),
                 corr0  = c(0.10, 0.08, 0.12, 0.09, 0.15, 0.07)),
    1.0, tolerance = tol
  )
})

# ---- guard test: correlations large enough to push adjusted rate toward 1 --

test_that("TR still sums to ~1 when large correlations trigger the clamp guard", {
  N <- 2
  # corr = 0.95 with rate01 = 0.07 would push adjusted rate to 0.07/0.05 = 1.4 > 1
  # The guard should clamp it and TR should still be close to 1
  expect_equal(
    tr_sum_check(bc = 0L,
                 rate10 = c(0.05, 0.08),
                 rate01 = c(0.07, 0.06),
                 corr1  = 0.0,
                 corr0  = 0.95),
    1.0, tolerance = 0.01   # wider tolerance: guard is an approximation at the boundary
  )
})
