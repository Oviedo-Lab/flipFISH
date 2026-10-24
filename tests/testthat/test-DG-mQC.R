
# File written by Claude Sonnet 5

# Tests for mQC_DG()/misread.qc.DG() -- the C++-backed DG/MCMCSA method, as distinct from the
# verbatim-preserved misread.qc.DG.legacy(). Three things are checked:
#   1. The core bivariate-normal / correlation-inversion math (bvn_cdf_uppertail_test(),
#      dg_find_sigma_RootBisection_test()), against R's own mvtnorm::pmvnorm and against exact
#      achieved-vs-target correlation recovery -- the actual bug misread.qc.DG() fixes relative to
#      misread.qc.DG.legacy() (see dg_sigma_formula()'s derivation in src/main.cpp).
#   2. A smoke test that misread.qc.DG() runs end-to-end and returns the expected shape, same
#      bundled sample data and output contract test-DG.R already checks for the legacy function.
#   3. A loose sanity cross-check against misread.qc.DG.legacy() on the same small sample (not an
#      exact-agreement check -- the two methods' internals now genuinely differ).

library(flipFISH)

tol <- 1e-8

# ---- 1. Core DG math --------------------------------------------------------------------------

test_that("bvn_cdf_uppertail matches R's mvtnorm::pmvnorm", {
  skip_if_not_installed("mvtnorm")
  hs   <- c(-1.0, -0.3, 0.2, 1.5)
  ks   <- c(-0.8, 0.1, 0.9, 2.0)
  rhos <- c(-0.8, -0.3, 0.0, 0.4, 0.9)
  for (h in hs) {
    for (k in ks) {
      for (rho in rhos) {
        sigma    <- matrix(c(1, rho, rho, 1), 2, 2)
        expected <- as.numeric(mvtnorm::pmvnorm(lower = c(h, k), upper = c(Inf, Inf), sigma = sigma))
        actual   <- bvn_cdf_uppertail_test(h, k, rho)
        expect_equal(actual, expected, tolerance = 1e-6)
      }
    }
  }
})

test_that("dg_find_sigma_RootBisection recovers the exact target correlation when feasible", {
  # Two bits with different (asymmetric) marginal rates, so threshold_i != threshold_j -- this is
  # exactly the case the single-shared-threshold historical/neuronsDG formula cannot handle.
  rate10 <- 0.03; rate01 <- 0.08
  t_i <- qnorm(rate10)        # bit i: true value 1, read=1 iff Z > t_i
  t_j <- qnorm(1 - rate01)    # bit j: true value 0, read=1 iff Z > t_j
  p_i <- 1 - pnorm(t_i)       # P(read_i = 1)
  p_j <- 1 - pnorm(t_j)       # P(read_j = 1)
  sd_i <- sqrt(p_i * (1 - p_i))
  sd_j <- sqrt(p_j * (1 - p_j))

  # Negative targets are feasible for this asymmetric pair (p_i near 1, p_j near 0); positive
  # targets beyond the Frechet-Hoeffding bound for these marginals are not, and should fall back to
  # rho = 0 (the achieved correlation then reflects independence, not the (infeasible) target).
  feasible_targets <- c(-0.5, -0.2, 0.0)
  for (rho_target in feasible_targets) {
    desired_cov   <- rho_target * sd_i * sd_j
    rho_gauss     <- dg_find_sigma_RootBisection_test(t_i, t_j, desired_cov)
    achieved_cov  <- bvn_cdf_uppertail_test(t_i, t_j, rho_gauss) - p_i * p_j
    achieved_corr <- achieved_cov / (sd_i * sd_j)
    expect_equal(achieved_corr, rho_target, tolerance = 1e-6)
  }

  infeasible_cov <- 0.6 * sd_i * sd_j
  expect_equal(dg_find_sigma_RootBisection_test(t_i, t_j, infeasible_cov), 0.0, tolerance = tol)
})

# ---- 2. misread.qc.DG() smoke test ------------------------------------------------------------

extdata_path <- system.file("extdata", "dg_legacy", package = "flipFISH")

load_small_STdata <- function(n_genes = 8, n_blanks = 3) {
  codebook <- read.csv(file.path(extdata_path, "codebook_obs.csv"), row.names = 1)
  genes    <- read.csv(file.path(extdata_path, "summary_stats_m1_genes.csv"))
  blanks   <- read.csv(file.path(extdata_path, "summary_stats_m1_blanks.csv"))
  blank_mask <- grepl("Blank", rownames(codebook), ignore.case = FALSE)
  rownames(genes)  <- rownames(codebook)[!blank_mask]
  rownames(blanks) <- rownames(codebook)[blank_mask]
  STdata <- rbind(genes[seq_len(n_genes), , drop = FALSE], blanks[seq_len(n_blanks), , drop = FALSE])
  list(STdata = STdata, codebook = codebook[rownames(STdata), , drop = FALSE])
}

test_that("misread.qc.DG runs end-to-end on a small subset of the bundled legacy data", {
  skip_if(extdata_path == "", "bundled DG legacy sample data not found")

  d <- load_small_STdata(n_genes = 8, n_blanks = 3)

  set.seed(1)
  res <- misread.qc.DG(
    d$STdata, d$codebook,
    resamples = 2, burnin = 5, FR_search = 5,
    n_tracker_updates = 2, n_forks = 1, N_cells = 100
  )

  expect_type(res, "list")
  expect_named(res, c("STdata", "fliprates", "erctc_plus", "msle", "PPV_by_step"))

  n_steps <- 2 + 5 + 5
  expect_equal(nrow(res$STdata), 11) # 8 genes + 3 blanks, echoing the (reordered) input
  expect_true(setequal(rownames(res$STdata), rownames(d$STdata)))
  n_bits <- ncol(d$codebook)
  expect_equal(length(res$fliprates), 2 * n_bits + n_bits * (n_bits - 1))
  expect_true(all(grepl("^rate10_bit|^rate01_bit|^corr_", names(res$fliprates))))
  expect_true(all(res$fliprates[grepl("^rate", names(res$fliprates))] >= 0 &
                  res$fliprates[grepl("^rate", names(res$fliprates))] <= 1))
  expect_equal(length(res$msle), n_steps)
  expect_equal(nrow(res$PPV_by_step), n_steps)
  expect_equal(ncol(res$PPV_by_step), 8)
  # check_data() sorts genes/blanks by decreasing observed rate before fitting (same convention
  # misread.qc() already uses), so column/row order need not match the caller's original row order --
  # just confirm it's the same set of names.
  expect_true(setequal(colnames(res$PPV_by_step), rownames(d$STdata)[seq_len(8)]))
  expect_named(res$erctc_plus, c("erc", "ecc", "ehc", "CR", "PPV", "PPV_lower", "PPV_upper"))
  expect_equal(nrow(res$erctc_plus), 11) # 8 genes + 3 blanks, same shape as misread.qc()'s own erctc_plus
  expect_true(setequal(rownames(res$erctc_plus), rownames(d$STdata)))

  expect_true(all(is.finite(res$msle)))
  expect_true(all(res$erctc_plus$PPV >= 0 & res$erctc_plus$PPV <= 1, na.rm = TRUE))
})

test_that("misread.qc.DG's feasibility penalty runs without error and stays finite", {
  skip_if(extdata_path == "", "bundled DG legacy sample data not found")
  d <- load_small_STdata(n_genes = 6, n_blanks = 3)

  set.seed(1)
  res <- misread.qc.DG(
    d$STdata, d$codebook,
    resamples = 2, burnin = 5, FR_search = 5,
    n_tracker_updates = 2, n_forks = 1, N_cells = 100,
    weights = list(prior = 0.1, dispersion = 0.1, feasibility = 0.05),
    max_fs_bc = 5
  )
  expect_true(all(is.finite(res$msle)))
})

# ---- 3. Second smoke test, different schedule/codebook subset ---------------------------------

# NOTE on why there is no automated misread.qc.DG() vs. misread.qc.DG.legacy() cross-check here:
# a live cross-check means running misread.qc.DG.legacy()'s full, unmodified MCMCSA random walk,
# and that walk's runtime is NOT reliably bounded even at a seemingly reasonable N_cells. During
# development, the same small codebook subset and schedule used elsewhere in this file ran
# misread.qc.DG.legacy() for over an hour with zero progress at N_cells=100 on one seed, while
# finishing in ~1 second on another -- Monte-Carlo-based simulated annealing's acceptance rate can
# collapse once the annealing temperature cools, and how often/how badly that happens is sensitive
# to the random draws, not just the nominal settings. misread.qc.DG() (this file's actual subject)
# does not show this on its own in the configurations tested above; but legacy is unmodified,
# verbatim-preserved code (see R/dg_legacy.R's own header), not something this change is meant to
# fix, and pinning a specific seed here would only be testing "does this one seed happen to be
# fast enough today" rather than anything real. A one-off manual comparison (both methods, same
# seed, same small codebook subset) was run during development and landed in the same broad PPV
# neighborhood for both; that comparison is not repeated here as an automated test.

test_that("misread.qc.DG runs end-to-end on a different small codebook subset/schedule", {
  skip_if(extdata_path == "", "bundled DG legacy sample data not found")
  d <- load_small_STdata(n_genes = 6, n_blanks = 3)

  set.seed(42)
  res <- misread.qc.DG(
    d$STdata, d$codebook,
    resamples = 2, burnin = 5, FR_search = 5,
    n_tracker_updates = 2, n_forks = 1, N_cells = 100
  )

  med <- median(res$erctc_plus$PPV, na.rm = TRUE)
  expect_true(is.finite(med) && med >= 0 && med <= 1)
  expect_true(all(is.finite(res$msle)))
})
