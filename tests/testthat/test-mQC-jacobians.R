
# File written by Claude Sonnet 5

# Tests for the single-point diagnostic exports mQC_ecc_jacobian_at(), mQC_ehc_jacobian_at(), and
# mQC_counts_at() (see their own comments in src/main.cpp, just above prepare_STdata_at()'s). These
# three were added so the analytic d(ecc)/d(theta) and d(ehc)/d(theta) Jacobians that
# expected_bc_counts() computes internally for mQC()'s L-BFGS fit and its delta-method PPV CI could
# be checked from R against finite differences, rather than only trusted via optimizer behavior --
# that intended use is exactly what's exercised below: mQC_counts_at()'s raw (erc, ecc, ehc) values
# at FR+-eps are differenced to build a numerical Jacobian, then compared to mQC_ecc_jacobian_at()/
# mQC_ehc_jacobian_at()'s analytic one.
#
# Until now, none of these three were ever actually called by a test (or by R/main.R or a vignette),
# so this finite-difference check had never been run. Doing so here surfaced a real bug: grep_idx()
# (src/main.cpp) returned the sentinel {-1} instead of an empty vector when a codebook had no
# "Blank"-labeled barcodes, and every caller of cb$blanks/cb$genes indexes into barcode-length
# vectors with it directly -- so a blank-free codebook caused an out-of-bounds write (confirmed with
# AddressSanitizer), whose effect depended on heap layout and so looked like nondeterministic noise
# in the Jacobians. Fixed by returning an empty vector instead; see git history on grep_idx(). The
# fixture below still includes two "Blank"-labeled rows precisely so this bug class stays covered.

library(flipFISH)

# ---- fixture: small synthetic codebook, tie-free under maxHam = 1 -----------------------------

# 4 barcodes, 6 bits, minimum pairwise Hamming distance 4 (so maxHam = 1 has no ambiguous/tied
# corrections to worry about) -- kept small so the O(n_params) finite-difference loop below is fast.
cb_mat <- matrix(c(
  0, 0, 1, 0, 1, 1,
  0, 1, 0, 1, 0, 1,
  1, 0, 0, 1, 1, 0,
  1, 1, 1, 0, 0, 0
), nrow = 4, byrow = TRUE)
rownames(cb_mat) <- c("gene1", "gene2", "Blank1", "Blank2")
storage.mode(cb_mat) <- "integer"

n_bits    <- ncol(cb_mat)
n_corr    <- n_bits * (n_bits - 1) / 2
maxHam    <- 1
n_forks   <- 1
max_flips <- n_bits

bc_counts <- cbind(
  rates    = c(500, 800, 50, 30),
  variance = c(5, 5, 5, 5),
  counts   = c(500, 800, 50, 30)
)
rownames(bc_counts) <- rownames(cb_mat)

# Moderate, nonzero flip rates/correlations -- order is [rate10 | rate01 | corr1 | corr0], matching
# pack_fr()'s layout.
set.seed(1)
FR0 <- c(
  rep(0.03, n_bits) + runif(n_bits, -0.005, 0.005),
  rep(0.07, n_bits) + runif(n_bits, -0.01,  0.01),
  runif(n_corr, -0.15, 0.15),
  runif(n_corr, -0.15, 0.15)
)
n_params <- length(FR0)

# Central finite difference of mQC_counts_at()'s requested column, evaluated at FR0.
fd_jacobian <- function(column, eps = 1e-6) {
  J <- matrix(0, nrow(bc_counts), n_params)
  for (k in seq_len(n_params)) {
    FRp <- FR0; FRp[k] <- FRp[k] + eps
    FRm <- FR0; FRm[k] <- FRm[k] - eps
    cp <- flipFISH:::mQC_counts_at(FRp, bc_counts, cb_mat, maxHam, n_forks, max_flips)
    cm <- flipFISH:::mQC_counts_at(FRm, bc_counts, cb_mat, maxHam, n_forks, max_flips)
    J[, k] <- (cp[, column] - cm[, column]) / (2 * eps)
  }
  J
}

# ---- mQC_counts_at() itself ---------------------------------------------------------------------

test_that("mQC_counts_at returns finite, non-negative erc/ecc/ehc with ehc <= ecc", {
  counts <- flipFISH:::mQC_counts_at(FR0, bc_counts, cb_mat, maxHam, n_forks, max_flips)
  expect_equal(colnames(counts), c("erc", "ecc", "ehc"))
  expect_equal(nrow(counts), nrow(bc_counts))
  expect_true(all(is.finite(counts)))
  expect_true(all(counts >= 0))
  expect_true(all(counts[, "ehc"] <= counts[, "ecc"] + 1e-8))
  # Blanks' true counts are zeroed by the misread-loss rescaling, so a blank can never itself be the
  # "hit" (true-barcode-is-BOI) term -- any nonzero ecc for a blank comes entirely from other
  # barcodes' misreads correcting to it.
  blank_rows <- grepl("Blank", rownames(bc_counts))
  expect_true(all(counts[blank_rows, "ehc"] == 0))
})

# ---- finite-difference checks -------------------------------------------------------------------

test_that("mQC_ecc_jacobian_at matches finite differences of mQC_counts_at's ecc column", {
  J_analytic <- flipFISH:::mQC_ecc_jacobian_at(FR0, bc_counts, cb_mat, maxHam, n_forks, max_flips)
  J_numeric  <- fd_jacobian("ecc")
  expect_equal(dim(J_analytic), c(nrow(bc_counts), n_params))
  expect_true(max(abs(J_analytic)) > 1) # sanity: fixture isn't accidentally flat/uninformative
  expect_equal(J_analytic, J_numeric, tolerance = 1e-5)
})

test_that("mQC_ehc_jacobian_at matches finite differences of mQC_counts_at's ehc column", {
  J_analytic <- flipFISH:::mQC_ehc_jacobian_at(FR0, bc_counts, cb_mat, maxHam, n_forks, max_flips)
  J_numeric  <- fd_jacobian("ehc")
  expect_equal(dim(J_analytic), c(nrow(bc_counts), n_params))
  expect_true(max(abs(J_analytic)) > 1)
  expect_equal(J_analytic, J_numeric, tolerance = 1e-5)
})
