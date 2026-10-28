# File written by Claude Sonnet 5

# Regression coverage for the boost::random draw sequences used throughout the C++ side (engine
# boost::random::mt19937, boost::random::normal_distribution, boost::random::bernoulli_distribution --
# see src/main.cpp's "RNG regression-test helpers" section). These are not statistical tests: the
# point is not whether the draws look normal/Bernoulli, but whether a *fixed seed* still produces
# the exact same numbers it did when these values were pinned. Package functions such as
# draw_random_fr() take an explicit seed for reproducibility; if a future change (reverting to
# std::, a Boost upgrade that swaps a distribution's internal sampling algorithm, or a refactor of
# call order/seeding) silently changes what a given seed produces, every downstream simulation
# result tied to that seed changes too, with no error anywhere else to catch it. These tests exist
# to make that change loud instead of silent.

tol <- 1e-12

test_that("mt19937 engine draw sequence for seed 12345 is unchanged", {
  expect_equal(
    flipFISH:::rng_engine_draw_test(12345L, 5L),
    c(3992670690, 3823185381, 1358822685, 561383553, 789925284)
  )
})

test_that("normal_distribution(0, 1) draw sequence for seed 12345 is unchanged", {
  expect_equal(
    flipFISH:::rng_normal_draw_test(12345L, 5L, 0.0, 1.0),
    c(-0.0967588508161653, 0.112274718607326, -0.410492815889905,
      -0.166521227996262, -1.80837181404269),
    tolerance = tol
  )
})

test_that("normal_distribution(2.5, 0.7) draw sequence for seed 42 is unchanged", {
  expect_equal(
    flipFISH:::rng_normal_draw_test(42L, 5L, 2.5, 0.7),
    c(2.052900378526, 1.91423431181635, 2.21960349339131, 1.89153781714527,
      1.8195687906043),
    tolerance = tol
  )
})

test_that("bernoulli_distribution(0.3) draw sequence for seed 7 is unchanged", {
  expect_equal(
    flipFISH:::rng_bernoulli_draw_test(7L, 10L, 0.3),
    c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)
  )
})

test_that("bernoulli_distribution(0.5) draw sequence for seed 999 is unchanged", {
  expect_equal(
    flipFISH:::rng_bernoulli_draw_test(999L, 10L, 0.5),
    c(0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L)
  )
})

test_that("normal_distribution draws via param_type reset (seed 2024) are unchanged", {
  # Mirrors draw_random_fr()'s reuse of one distribution object across different parameterizations
  # (its normCorr / sdcorr_pairs path) instead of constructing a fresh distribution per draw.
  expect_equal(
    flipFISH:::rng_normal_paramtype_draw_test(2024L, 0.0, 1.0, 5.0, 2.0),
    c(-0.084637885913434, 1.79729508341745),
    tolerance = tol
  )
})

# gamma_distribution/poisson_distribution: added when dg_draw_true_counts()'s per-cell transcript
# draws moved from R::rgamma()/R::rpois() to boost::random, so mQC_DG()'s hot loop no longer calls
# into R at all (see dg_draw_true_counts(), and the `seed` argument added to mQC_DG()/
# misread.qc.DG()).

test_that("gamma_distribution(2.0, 1.5) draw sequence for seed 12345 is unchanged", {
  expect_equal(
    flipFISH:::rng_gamma_draw_test(12345L, 5L, 2.0, 1.5),
    c(0.915970170923289, 5.49268773199246, 3.1940126259417,
      3.44508430664966, 1.14114138563479),
    tolerance = tol
  )
})

test_that("gamma_distribution(5.0, 0.5) draw sequence for seed 42 is unchanged", {
  expect_equal(
    flipFISH:::rng_gamma_draw_test(42L, 5L, 5.0, 0.5),
    c(1.76587267434506, 2.8003654916141, 2.80021950625805,
      2.27679153234628, 1.32938534881185),
    tolerance = tol
  )
})

test_that("poisson_distribution(3.2) draw sequence for seed 7 is unchanged", {
  expect_equal(
    flipFISH:::rng_poisson_draw_test(7L, 10L, 3.2),
    c(1L, 2L, 4L, 2L, 3L, 7L, 4L, 3L, 7L, 2L)
  )
})

test_that("poisson_distribution(0.5) draw sequence for seed 999 is unchanged", {
  expect_equal(
    flipFISH:::rng_poisson_draw_test(999L, 10L, 0.5),
    c(1L, 0L, 0L, 3L, 0L, 0L, 1L, 0L, 0L, 0L)
  )
})

# dg_ppois_lower(): this file's native (boost::math, no R callback) replacement for the single
# R::ppois(x, lambda, 1, 0) call dg_draw_true_counts() used to make per MCMCSA step -- a statistical
# correctness check against R's own ppois(), not a pinned-sequence test (no randomness involved).

test_that("dg_ppois_lower matches R's ppois(x, lambda, lower.tail = TRUE) over a grid", {
  xs      <- c(0, 0.5, 1, 2.5, 3, 7, 10.999, 15)
  lambdas <- c(0.5, 1, 3.2, 10, 25)
  for (lambda in lambdas) {
    for (x in xs) {
      expect_equal(flipFISH:::dg_ppois_lower_test(x, lambda), ppois(x, lambda), tolerance = 1e-9)
    }
  }
})

test_that("dg_ppois_lower reproduces R's ppois() edge cases", {
  expect_equal(flipFISH:::dg_ppois_lower_test(5, 0), 1) # lambda == 0: Poisson(0) is a point mass at 0
  expect_equal(flipFISH:::dg_ppois_lower_test(-1, 3), 0) # x < 0
})
