args <- commandArgs(trailingOnly = TRUE)
mouse_file <- args[1]   # e.g. "m1"
maxHam_val <- as.integer(args[2])
out_path   <- args[3]

library(flipFISH)
base <- "old_vs_new_ppv_comparison/input_data"
codebook <- read.csv(file.path(base, "codebook_obs.csv"), row.names = 1)
g <- read.csv(file.path(base, "summary_stats", paste0("summary_stats_", mouse_file, "_genes.csv")))
b <- read.csv(file.path(base, "summary_stats", paste0("summary_stats_", mouse_file, "_blanks.csv")))
rownames(g) <- rownames(codebook)[1:721]
rownames(b) <- rownames(codebook)[722:819]
STdata <- rbind(g, b)

t0 <- Sys.time()
qc <- misread.qc(
  STdata, codebook,
  n_forks   = 20,
  max_flips = 6,
  maxeval   = 500,
  maxHam    = maxHam_val,
  weights   = list(blank = 10, prior = 0.5, dispersion = 0.5, erc = 0)
)
t1 <- Sys.time()
cat("\nElapsed:", as.numeric(difftime(t1, t0, units = "secs")), "sec\n")

saveRDS(qc, out_path)
cat("Saved to", out_path, "\n")
