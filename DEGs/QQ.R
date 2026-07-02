#!/usr/bin/env Rscript
################################################################################
# Create QQ plots comparing real data to permutation null distribution
################################################################################

library(data.table)
library(ggplot2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript create_permutation_QQ_plots.R <timepoint> <n_perms>\n  timepoint: 6 or 24\n  n_perms: number of permutations")
}

tp      <- args[1]
n_perms <- as.integer(args[2])

# Configuration
base_dir   <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates"
deg_dir    <- file.path(base_dir, "DEG")
perm_dir   <- file.path(base_dir, "DEG", "permutations")
output_dir <- file.path(deg_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
today <- format(Sys.time(), "%m%d%Y")

cat("\n=== Creating QQ plots ===\n")
cat("Timepoint: T", tp, "\n")

################################################################################
# Load real data
################################################################################

cat("\n=== Loading real data ===\n")

real_file <- tail(sort(list.files(deg_dir,
  pattern = paste0(".*_MBP_T", tp, "\\.txt$"),
  full.names = TRUE)), 1)

if (length(real_file) == 0 || !file.exists(real_file)) {
  stop("Real data file not found in: ", deg_dir)
}
cat("Real data file:", basename(real_file), "\n")

real_data <- fread(real_file)
if (colnames(real_data)[1] %in% c("V1", "", "Gene", "gene", "gene_id")) {
  setnames(real_data, 1, "ensgene")
}
cat("Loaded real data:", nrow(real_data), "genes\n")

real_pvals <- real_data$P.Value[!is.na(real_data$P.Value)]
cat("Non-NA p-values:", length(real_pvals), "\n")

################################################################################
# Load permutation files
################################################################################

cat("\n=== Loading permutation data ===\n")

perm_files <- sort(list.files(perm_dir,
  pattern = paste0(".*_MBP_T", tp, "_PERM[0-9]+\\.txt$"),
  full.names = TRUE))

if (length(perm_files) == 0) stop("No permutation files found in: ", perm_dir)
cat("Found", length(perm_files), "permutation files\n")

################################################################################
# QQ data helper
################################################################################

get_qq_data <- function(pvals) {
  n        <- length(pvals)
  observed <- -log10(sort(pvals))
  expected <- -log10((1:n) / (n + 1))
  data.frame(expected = expected, observed = observed)
}

################################################################################
# Load all permutation p-values
################################################################################

cat("\n=== Processing permutation p-values ===\n")

perm_qq_list <- vector("list", length(perm_files))

for (i in seq_along(perm_files)) {
  if (i %% 20 == 0) cat("Processing permutation", i, "/", length(perm_files), "\n")
  perm_data        <- fread(perm_files[i])
  perm_pvals       <- perm_data$P.Value[!is.na(perm_data$P.Value)]
  qq_data          <- get_qq_data(perm_pvals)
  qq_data$perm     <- i
  perm_qq_list[[i]] <- qq_data
}

perm_qq_all <- do.call(rbind, perm_qq_list)
cat("Total permutation points:", nrow(perm_qq_all), "\n")

################################################################################
# Create QQ plot
################################################################################

cat("\n=== Creating QQ plot ===\n")

real_qq <- get_qq_data(real_pvals)

p <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "black", linewidth = 0.8) +
  geom_line(data = perm_qq_all,
            aes(x = expected, y = observed, group = perm),
            color = "gray70", alpha = 0.3, linewidth = 0.5) +
  geom_line(data = real_qq,
            aes(x = expected, y = observed),
            color = "red", linewidth = 1.5) +
  labs(
    title    = paste("MBP vs EtOH T", tp, "- QQ Plot"),
    subtitle = paste("Real data (red) vs", length(perm_files), "permutations (gray)"),
    x        = "Expected -log10(p)",
    y        = "Observed -log10(p)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    axis.ticks       = element_line(color = "black"),
    plot.title       = element_text(face = "bold", size = 16, color = "black"),
    plot.subtitle    = element_text(size = 12, color = "black")
  )

output_file <- file.path(output_dir,
                         paste0(today, "_MBP_T", tp, "_perm_QQ_plot.png"))
ggsave(output_file, p, width = 10, height = 10, dpi = 300)
cat("Saved:", output_file, "\n")

################################################################################
# Genomic inflation factor (lambda)
################################################################################

cat("\n=== Genomic Inflation Factor ===\n")

calc_lambda <- function(pvals) {
  chisq <- qchisq(1 - pvals, 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
}

lambda_real <- calc_lambda(real_pvals)
cat("Real data lambda:", round(lambda_real, 3), "\n")

lambda_perms <- sapply(perm_files, function(f) {
  perm_data  <- fread(f)
  perm_pvals <- perm_data$P.Value[!is.na(perm_data$P.Value)]
  calc_lambda(perm_pvals)
})

cat("Permutation lambda:\n")
cat("  Mean:",   round(mean(lambda_perms),   3), "\n")
cat("  Median:", round(median(lambda_perms), 3), "\n")
cat("  SD:",     round(sd(lambda_perms),     3), "\n")
cat("  Range:",  round(min(lambda_perms),    3), "-", round(max(lambda_perms), 3), "\n")

lambda_summary <- data.frame(
  type   = c("real", rep("permutation", length(lambda_perms))),
  lambda = c(lambda_real, lambda_perms)
)

lambda_file <- file.path(output_dir,
                         paste0(today, "_MBP_T", tp, "_lambda_summary.txt"))
write.table(lambda_summary, lambda_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Saved lambda summary:", lambda_file, "\n")

cat("\n=== COMPLETE ===\n")
