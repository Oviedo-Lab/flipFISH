
# File written by Claude Sonnet 5

# Smoke test for misread.qc.DG.legacy(), the verbatim-ported legacy dichotomized-Gaussian (DG)
# spot-simulation method (see R/dg_legacy.R, src/dg_legacy.cpp). Uses a small subset of the bundled m1
# sample data (inst/extdata/dg_legacy/) -- a tiny codebook/cell count just to confirm the ported pipeline
# runs end-to-end and returns results of the expected shape; it is not a correctness/recovery check.

library(flipFISH)

extdata_path <- system.file("extdata", "dg_legacy", package = "flipFISH")

test_that("misread.qc.DG.legacy runs end-to-end on a small subset of the bundled legacy data", {
  skip_if(extdata_path == "", "bundled DG legacy sample data not found")

  codebook <- read.csv(file.path(extdata_path, "codebook_obs.csv"), row.names = 1)
  genes    <- read.csv(file.path(extdata_path, "summary_stats_m1_genes.csv"))
  blanks   <- read.csv(file.path(extdata_path, "summary_stats_m1_blanks.csv"))

  # The bundled summary-stats CSVs have no row names: their rows are positionally aligned with the
  # codebook's rows (genes first, then blanks) -- see the "Reorder codebook" note in misread.qc.DG().
  blank_mask <- grepl("Blank", rownames(codebook), ignore.case = FALSE)
  rownames(genes)  <- rownames(codebook)[!blank_mask]
  rownames(blanks) <- rownames(codebook)[blank_mask]

  # Subset down to a small, fast test problem (a handful of genes/blanks, not the full ~819-barcode
  # codebook) so this test runs in a couple of seconds rather than several minutes.
  n_genes  <- 15
  n_blanks <- 5
  STdata   <- rbind(genes[seq_len(n_genes), , drop = FALSE], blanks[seq_len(n_blanks), , drop = FALSE])
  codebook_small <- codebook[rownames(STdata), , drop = FALSE]

  set.seed(1)
  res <- misread.qc.DG.legacy(
    STdata, codebook_small,
    resamples = 2, burnin = 5, FR_search = 5,
    n_tracker_updates = 2, n_forks = 1, N_cells = 20
  )

  expect_type(res, "list")
  expect_named(
    res,
    c("sim_summaries", "PPV", "PPV_summary", "Counts_sorted_sim", "Counts_sorted_obs", "PPV_plot", "sorted_count_plot")
  )

  n_steps <- 2 + 5 + 5
  expect_equal(nrow(res$sim_summaries), n_steps)
  expect_equal(nrow(res$PPV), n_steps)
  expect_equal(ncol(res$PPV), n_genes)
  expect_equal(nrow(res$PPV_summary), n_genes)
  expect_equal(nrow(res$Counts_sorted_sim), n_steps)
  expect_equal(ncol(res$Counts_sorted_sim), n_genes + n_blanks)

  # sim_summaries legitimately contains NA in some columns at this small a sample size (e.g.
  # PPV_mean_filtered/PPV_sd_filtered when the "filtered" gene group has 0-1 members; flip_rate_10/
  # flip_rate_01, the back-calculated-from-spots columns, are NA except on the few steps around
  # FR_search where return_flip_rates is toggled on -- see .DG_run_simulation()) -- check the columns
  # that should always be defined instead of requiring the whole data frame to be finite.
  expect_true(all(is.finite(res$sim_summaries$CR_PPV_slope)))
  expect_true(all(is.finite(res$sim_summaries$sim_mse)))
  expect_true(all(res$PPV_summary$median >= 0 & res$PPV_summary$median <= 1, na.rm = TRUE))

  expect_s3_class(res$PPV_plot, "ggplot")
  expect_s3_class(res$sorted_count_plot, "ggplot")
})
