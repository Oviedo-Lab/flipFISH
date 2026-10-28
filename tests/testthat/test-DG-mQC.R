
# File written by Claude Sonnet 5

# Tests for mQC_DG()/misread.qc.DG() -- the C++-backed DG/MCMCSA method. Three things are checked:
#   1. The core bivariate-normal / correlation-inversion math (bvn_cdf_uppertail_test(),
#      dg_find_sigma_RootBisection_test()), against R's own mvtnorm::pmvnorm and against exact
#      achieved-vs-target correlation recovery (see dg_sigma_formula()'s derivation in src/main.cpp).
#   2. Smoke tests that misread.qc.DG() runs end-to-end and returns the expected shape, on a small
#      subset of the bundled sample data (inst/extdata/sample_m1/).
#   3. That max_calls_per_step actually bounds worst-case runtime regardless of seed (see its own
#      test below): simulated annealing's acceptance rate can collapse for some (seed, data,
#      schedule) combinations, needing many proposal evaluations to advance a single step even when
#      each evaluation itself is cheap -- not just a hypothetical, this is exactly what made these
#      smoke tests intermittently take minutes before max_calls_per_step was added (see git log).
#
# misread.qc.DG()'s own randomness (proposals, stochastic pruning, acceptance, per-cell transcript
# draws) is driven entirely by its own `seed` argument via a dedicated boost::random::mt19937 on the
# C++ side, not by R's RNG -- so set.seed() has no effect on it, and these tests don't call it.

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
        actual   <- flipFISH:::bvn_cdf_uppertail_test(h, k, rho)
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
    rho_gauss     <- flipFISH:::dg_find_sigma_RootBisection_test(t_i, t_j, desired_cov)
    achieved_cov  <- flipFISH:::bvn_cdf_uppertail_test(t_i, t_j, rho_gauss) - p_i * p_j
    achieved_corr <- achieved_cov / (sd_i * sd_j)
    expect_equal(achieved_corr, rho_target, tolerance = 1e-6)
  }

  infeasible_cov <- 0.6 * sd_i * sd_j
  expect_equal(flipFISH:::dg_find_sigma_RootBisection_test(t_i, t_j, infeasible_cov), 0.0, tolerance = tol)
})

# ---- 2. misread.qc.DG() smoke test ------------------------------------------------------------

extdata_path <- system.file("extdata", "sample_m1", package = "flipFISH")

load_small_STdata <- function(gene_idx, n_blanks = 3) {
  codebook <- read.csv(file.path(extdata_path, "codebook_obs.csv"), row.names = 1)
  genes    <- read.csv(file.path(extdata_path, "summary_stats_m1_genes.csv"))
  blanks   <- read.csv(file.path(extdata_path, "summary_stats_m1_blanks.csv"))
  blank_mask <- grepl("Blank", rownames(codebook), ignore.case = FALSE)
  rownames(genes)  <- rownames(codebook)[!blank_mask]
  rownames(blanks) <- rownames(codebook)[blank_mask]
  STdata <- rbind(genes[gene_idx, , drop = FALSE], blanks[seq_len(n_blanks), , drop = FALSE])
  list(STdata = STdata, codebook = codebook[rownames(STdata), , drop = FALSE])
}

# gene_idx picks the package's own lowest-expression-rate genes (summary_stats_m1_genes.csv's
# "rates" column), not just the first few by row order: misread.qc.DG()'s MCMCSA per-call cost
# scales with how many transcripts get simulated per gene (proportional to each gene's rate, not to
# N_cells alone -- see dg_draw_true_counts()), and a shape/smoke test doesn't need realistic
# expression levels. This keeps each MCMCSA proposal evaluation these smoke tests make cheap. The
# schedule (resamples/burnin/FR_search) is also kept minimal for the same reason, and
# max_calls_per_step caps worst-case runtime outright: simulated annealing's acceptance rate can
# collapse for some (seed, data, schedule) combinations regardless of how cheap each proposal is
# (requiring many of them to advance a single step -- see misread.qc.DG()'s own doc on this), so
# cheap proposals alone don't guarantee a fast test, only a bounded number of them does.

test_that("misread.qc.DG runs end-to-end on a small subset of the bundled legacy data", {
  skip_if(extdata_path == "", "bundled sample_m1 data not found")

  d <- load_small_STdata(gene_idx = c(5, 8, 4), n_blanks = 3) # rates ~0.08/0.16/0.21

  res <- misread.qc.DG(
    d$STdata, d$codebook,
    resamples = 1, burnin = 2, FR_search = 2,
    n_tracker_updates = 0, n_forks = 1, N_cells = 100,
    max_calls_per_step = 10
  )

  expect_type(res, "list")
  expect_named(res, c("STdata", "fliprates", "erctc_plus", "msle", "PPV_by_step"))

  n_steps <- 1 + 2 + 2
  expect_equal(nrow(res$STdata), 6) # 3 genes + 3 blanks, echoing the (reordered) input
  expect_true(setequal(rownames(res$STdata), rownames(d$STdata)))
  n_bits <- ncol(d$codebook)
  expect_equal(length(res$fliprates), 2 * n_bits + n_bits * (n_bits - 1))
  expect_true(all(grepl("^rate10_bit|^rate01_bit|^corr_", names(res$fliprates))))
  expect_true(all(res$fliprates[grepl("^rate", names(res$fliprates))] >= 0 &
                  res$fliprates[grepl("^rate", names(res$fliprates))] <= 1))
  expect_equal(length(res$msle), n_steps)
  expect_equal(nrow(res$PPV_by_step), n_steps)
  expect_equal(ncol(res$PPV_by_step), 3)
  # check_data() sorts genes/blanks by decreasing observed rate before fitting (same convention
  # misread.qc() already uses), so column/row order need not match the caller's original row order --
  # just confirm it's the same set of names.
  expect_true(setequal(colnames(res$PPV_by_step), rownames(d$STdata)[seq_len(3)]))
  expect_named(res$erctc_plus, c("erc", "ecc", "ehc", "CR", "PPV", "PPV_lower", "PPV_upper"))
  expect_equal(nrow(res$erctc_plus), 6) # 3 genes + 3 blanks, same shape as misread.qc()'s own erctc_plus
  expect_true(setequal(rownames(res$erctc_plus), rownames(d$STdata)))

  expect_true(all(is.finite(res$msle)))
  expect_true(all(res$erctc_plus$PPV >= 0 & res$erctc_plus$PPV <= 1, na.rm = TRUE))
})

test_that("misread.qc.DG's feasibility penalty runs without error and stays finite", {
  skip_if(extdata_path == "", "bundled sample_m1 data not found")
  d <- load_small_STdata(gene_idx = c(5, 8), n_blanks = 3) # rates ~0.08/0.16

  res <- misread.qc.DG(
    d$STdata, d$codebook,
    resamples = 1, burnin = 2, FR_search = 2,
    n_tracker_updates = 0, n_forks = 1, N_cells = 100,
    weights = list(prior = 0.1, dispersion = 0.1, feasibility = 0.05),
    max_fs_bc = 1, max_calls_per_step = 10
  )
  expect_true(all(is.finite(res$msle)))
})

# ---- 3. max_calls_per_step bounds worst-case runtime regardless of seed -----------------------

test_that("max_calls_per_step bounds misread.qc.DG's runtime across seeds that would otherwise stall", {
  skip_if(extdata_path == "", "bundled sample_m1 data not found")
  # Same feasibility-penalty configuration test 2 uses (max_fs_bc = 1 barcode is already the
  # cheapest possible per-call cost), but with a longer schedule and max_calls_per_step = 1 -- which
  # forces every single step to accept its first proposal, so each of these runs needs exactly
  # n_steps = 15 total MCMCSA calls no matter how the annealing's acceptance rate behaves for a
  # given seed. Several seeds known to need many dozens of calls without a cap (see this file's git
  # history) are checked here specifically, to confirm the cap -- not seed luck -- is what's
  # keeping this fast. "Seed" here is misread.qc.DG()'s own `seed` argument (its C++-side
  # boost::random::mt19937), not R's set.seed() -- see this file's top-of-file note.
  d <- load_small_STdata(gene_idx = c(5, 8), n_blanks = 3)
  n_steps <- 5 + 5 + 5

  for (seed in c(1, 2, 6, 7, 9, 10)) {
    t0  <- Sys.time()
    res <- misread.qc.DG(
      d$STdata, d$codebook,
      resamples = 5, burnin = 5, FR_search = 5,
      n_tracker_updates = 0, n_forks = 1, N_cells = 100,
      weights = list(prior = 0.1, dispersion = 0.1, feasibility = 0.05),
      max_fs_bc = 1, max_calls_per_step = 1,
      seed = seed
    )
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    expect_lt(elapsed, 10) # generous margin; uncapped, these seeds take tens of seconds or more
    expect_equal(length(res$msle), n_steps)
    expect_true(all(is.finite(res$msle)))
  }
})

# ---- 4. Second smoke test, different schedule/codebook subset ---------------------------------

test_that("misread.qc.DG runs end-to-end on a different small codebook subset/schedule", {
  skip_if(extdata_path == "", "bundled sample_m1 data not found")
  d <- load_small_STdata(gene_idx = c(8, 4), n_blanks = 2) # rates ~0.16/0.21

  res <- misread.qc.DG(
    d$STdata, d$codebook,
    resamples = 1, burnin = 2, FR_search = 2,
    n_tracker_updates = 0, n_forks = 1, N_cells = 100,
    max_calls_per_step = 10
  )

  med <- median(res$erctc_plus$PPV, na.rm = TRUE)
  expect_true(is.finite(med) && med >= 0 && med <= 1)
  expect_true(all(is.finite(res$msle)))
})
