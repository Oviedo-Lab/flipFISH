## Compare misread.qc() results at maxHam=1 vs. maxHam=3 (both run fresh on current
## codebase, ACx-only real data, n_forks=20), and sanity-check the fresh maxHam=3 run
## against the previously-saved maxHam=3 results (old_vs_new_ppv_comparison/qc_outputs/),
## which were produced by an earlier package commit (3861d97) before a large src/main.cpp
## cleanup (9067c2f, ec4abcc).

setwd("/Users/michaelbarkasi/projects_Oviedo_lab/R_packages/flipFISH")

mice      <- c("m1", "m2", "m3", "m4", "m8")
mouse_num <- c(m1 = 1, m2 = 2, m3 = 3, m4 = 4, m8 = 5)  # old-run numbering, per README

new_dir <- "old_vs_new_ppv_comparison/rerun_maxHam1_vs_maxHam3/qc_outputs"
old_dir <- "old_vs_new_ppv_comparison/qc_outputs"

rows <- list()
for (m in mice) {
  qc3 <- readRDS(file.path(new_dir, paste0("qc_real_", m, "_maxHam3.rds")))
  qc1 <- readRDS(file.path(new_dir, paste0("qc_real_", m, "_maxHam1.rds")))
  old_rds <- file.path(old_dir, paste0("qc_real_mouse", mouse_num[m], ".rds"))
  qc_old  <- if (file.exists(old_rds)) readRDS(old_rds) else NULL

  genes <- rownames(qc3$erctc_plus)
  stopifnot(identical(genes, rownames(qc1$erctc_plus)))

  df <- data.frame(
    mouse        = mouse_num[m],
    mouse_label  = m,
    gene         = genes,
    is_blank     = grepl("Blank", genes),
    PPV_maxHam3  = qc3$erctc_plus$PPV,
    PPV_maxHam1  = qc1$erctc_plus$PPV,
    ecc_maxHam3  = qc3$erctc_plus$ecc,
    ecc_maxHam1  = qc1$erctc_plus$ecc,
    erc_maxHam3  = qc3$erctc_plus$erc,
    erc_maxHam1  = qc1$erctc_plus$erc,
    CR_maxHam3   = qc3$erctc_plus$CR,
    CR_maxHam1   = qc1$erctc_plus$CR
  )

  if (!is.null(qc_old)) {
    # Old results are in original fit order; align by gene name
    old_tab <- qc_old$erctc_plus
    df$PPV_old_maxHam3 <- old_tab[genes, "PPV"]
  } else {
    df$PPV_old_maxHam3 <- NA_real_
  }

  rows[[m]] <- df
}
full <- do.call(rbind, rows)
write.csv(full, "old_vs_new_ppv_comparison/rerun_maxHam1_vs_maxHam3/gene_level_comparison.csv", row.names = FALSE)

cat("================================================================\n")
cat("SANITY CHECK: fresh maxHam=3 (current code) vs. saved maxHam=3 (commit 3861d97)\n")
cat("================================================================\n")
if (all(!is.na(full$PPV_old_maxHam3))) {
  r_sanity  <- cor(full$PPV_maxHam3, full$PPV_old_maxHam3, method = "pearson")
  rs_sanity <- cor(full$PPV_maxHam3, full$PPV_old_maxHam3, method = "spearman")
  max_abs_diff <- max(abs(full$PPV_maxHam3 - full$PPV_old_maxHam3))
  med_abs_diff <- median(abs(full$PPV_maxHam3 - full$PPV_old_maxHam3))
  cat(sprintf("Pooled (n=%d): Pearson r = %.6f, Spearman r = %.6f\n", nrow(full), r_sanity, rs_sanity))
  cat(sprintf("Median |delta PPV| = %.6f, Max |delta PPV| = %.6f\n", med_abs_diff, max_abs_diff))
  for (m in mice) {
    sub <- full[full$mouse_label == m, ]
    cat(sprintf("  mouse %s: Pearson r = %.6f, median|delta| = %.6f, max|delta| = %.6f\n",
                m, cor(sub$PPV_maxHam3, sub$PPV_old_maxHam3), median(abs(sub$PPV_maxHam3 - sub$PPV_old_maxHam3)),
                max(abs(sub$PPV_maxHam3 - sub$PPV_old_maxHam3))))
  }
} else {
  cat("Old saved results missing for at least one mouse; skipping sanity check.\n")
}

cat("\n================================================================\n")
cat("MAIN COMPARISON: maxHam=1 vs. maxHam=3 (both fresh, current code)\n")
cat("================================================================\n")
r_main  <- cor(full$PPV_maxHam1, full$PPV_maxHam3, method = "pearson")
rs_main <- cor(full$PPV_maxHam1, full$PPV_maxHam3, method = "spearman")
cat(sprintf("Pooled (n=%d gene-mouse pairs): Pearson r = %.4f, Spearman r = %.4f\n", nrow(full), r_main, rs_main))

summary_rows <- list()
for (m in mice) {
  sub   <- full[full$mouse_label == m & !full$is_blank, ]  # genes only for "good gene" comparison
  good3 <- sub$PPV_maxHam3 >= 0.8
  good1 <- sub$PPV_maxHam1 >= 0.8
  both    <- sum(good3 & good1)
  only3   <- sum(good3 & !good1)
  only1   <- sum(!good3 & good1)
  jacc    <- both / (both + only3 + only1)
  summary_rows[[m]] <- data.frame(
    mouse = mouse_num[m], mouse_label = m, n_genes = nrow(sub),
    n_good_maxHam3 = sum(good3), n_good_maxHam1 = sum(good1),
    both_good = both, maxHam3_only = only3, maxHam1_only = only1,
    jaccard = jacc,
    pearson_r = cor(sub$PPV_maxHam1, sub$PPV_maxHam3),
    spearman_r = cor(sub$PPV_maxHam1, sub$PPV_maxHam3, method = "spearman"),
    median_ecc_ratio = median(sub$ecc_maxHam1 / sub$ecc_maxHam3)
  )
}
summary_df <- do.call(rbind, summary_rows)
print(summary_df, row.names = FALSE)
write.csv(summary_df, "old_vs_new_ppv_comparison/rerun_maxHam1_vs_maxHam3/summary_by_mouse.csv", row.names = FALSE)

cat("\nOverall good-gene counts (genes only, pooled):\n")
genes_only <- full[!full$is_blank, ]
cat(sprintf("  maxHam=3 good: %d / %d\n", sum(genes_only$PPV_maxHam3 >= 0.8), nrow(genes_only)))
cat(sprintf("  maxHam=1 good: %d / %d\n", sum(genes_only$PPV_maxHam1 >= 0.8), nrow(genes_only)))
cat(sprintf("  both good: %d\n", sum(genes_only$PPV_maxHam3 >= 0.8 & genes_only$PPV_maxHam1 >= 0.8)))
cat(sprintf("  maxHam3-only good: %d\n", sum(genes_only$PPV_maxHam3 >= 0.8 & genes_only$PPV_maxHam1 < 0.8)))
cat(sprintf("  maxHam1-only good: %d\n", sum(genes_only$PPV_maxHam3 < 0.8 & genes_only$PPV_maxHam1 >= 0.8)))

## Scatter plot
suppressMessages(library(ggplot2))
p <- ggplot(full, aes(x = PPV_maxHam3, y = PPV_maxHam1, color = is_blank)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "orange"), name = "Blank") +
  labs(x = "PPV (maxHam = 3)", y = "PPV (maxHam = 1)",
       title = "misread.qc() PPV: maxHam=1 vs. maxHam=3 (ACx-only, 5 mice, pooled)") +
  theme_minimal()
ggsave("old_vs_new_ppv_comparison/rerun_maxHam1_vs_maxHam3/maxHam1_vs_maxHam3_scatter.png",
       p, width = 6, height = 5, dpi = 150)

cat("\nWrote gene_level_comparison.csv, summary_by_mouse.csv, maxHam1_vs_maxHam3_scatter.png\n")
