set.seed(8214)
t0 <- Sys.time()
bm_erc <- sim.benchmark(
  m1_summary_stats,
  codebook,
  n_sims                           = 5,
  n_forks                          = 20,
  max_flips                        = 6,
  report_freq                      = 20,
  maxeval                          = 200,
  max_correctable_Hamming_distance = NULL,
  blank_weight                     = 10,
  prior_weight                     = 7,
  erc_weight                       = 1,    # or 3, per the "stronger emphasis" option
  n_restarts                       = 1     # or 3, per the "stronger emphasis" option
)
cat("\nElapsed (erc + multi-start, full scale):", round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min\n")

N_bits <- 28
l      <- length(bm_erc$fr_stip)   # same length for bm_reg2

build_metrics <- function(bm) {
  n_sims <- nrow(bm$ecc_est)
  
  # Bit-flip correlation terms: stipulated vs. estimated (all sims, all barcode pairs)
  cor_stip <- rep(bm$fr_stip[(2*N_bits + 1):l], each = n_sims)
  cor_est  <- c(bm$fr__est[, (2*N_bits + 1):l])
  cor_stip[is.na(cor_stip)] <- 0.0
  cor_est[is.na(cor_est)]   <- 0.0
  
  # PPV: stipulated (expected) vs. estimated, averaged across sims per barcode
  PPV_exp <- bm$PPV_exp
  PPV_est <- colMeans(bm$PPV_est, na.rm = TRUE)
  
  # Low-true-PPV band: barcodes with stipulated PPV in (0, 0.4)
  df_low <- subset(data.frame(PPV_exp = PPV_exp, PPV_est = PPV_est),
                   PPV_exp > 0 & PPV_exp < 0.4)
  
  list(
    cor_recovery       = cor(cor_stip, cor_est),
    mean_abs_corr_est  = mean(abs(cor_est)),
    corr_est_range     = range(cor_est),
    overall_ppv_cor    = cor(PPV_exp, PPV_est, use = "complete.obs"),
    low_ppv_band       = df_low,
    mean_low_ppv_est   = mean(df_low$PPV_est, na.rm = TRUE),
    mean_low_ppv_true  = mean(df_low$PPV_exp, na.rm = TRUE)
  )
}

m_base <- build_metrics(bm_base2)
m_reg  <- build_metrics(bm_reg2)
m_erc  <- build_metrics(bm_erc)

bias_base <- m_base$mean_low_ppv_est - m_base$mean_low_ppv_true
bias_reg  <- m_reg$mean_low_ppv_est  - m_reg$mean_low_ppv_true
bias_erc  <- m_erc$mean_low_ppv_est  - m_erc$mean_low_ppv_true

comparison_table <- data.frame(
  metric = c(
    "cor_recovery (bit-flip correlation terms)",
    "mean |corr estimate|",
    "overall cor(PPV_exp, PPV_est)",
    "mean PPV_est, barcodes with true PPV in (0.13-0.4)",
    "absolute bias in that low-PPV band"
  ),
  baseline = c(
    m_base$cor_recovery,
    m_base$mean_abs_corr_est,
    m_base$overall_ppv_cor,
    m_base$mean_low_ppv_est,
    bias_base
  ),
  regularized = c(
    m_reg$cor_recovery,
    m_reg$mean_abs_corr_est,
    m_reg$overall_ppv_cor,
    m_reg$mean_low_ppv_est,
    bias_reg
  ),
  erc = c(
    m_erc$cor_recovery,
    m_erc$mean_abs_corr_est,
    m_erc$overall_ppv_cor,
    m_erc$mean_low_ppv_est,
    bias_erc
  ),
  true_value = c(
    1.0,   # ideal
    mean(abs(bm_reg2$fr_stip[(2*N_bits+1):l])),  # ~0.125 mean |stipulated corr| computed separately, see note below
    1.0,   # ideal
    m_reg$mean_low_ppv_true,
    0.0    # ideal
  )
)

comparison_table