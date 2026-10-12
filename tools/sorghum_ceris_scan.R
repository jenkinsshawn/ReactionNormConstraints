#!/usr/bin/env Rscript
# =============================================================================
# sorghum_ceris_scan.R  (ReactionNormConstraints v10.0.0)
#
# Reimplements the CERIS search of Li et al. (2018) on their sorghum data and
# reports, for every candidate index, its correlation with the environment
# means (rho) and the largest possible loading on environmental variation
# unrelated to the environment means, tau_max = sqrt(1 - rho^2).
#
# Data: the public CERIS-JGRA workshop repository,
#   https://github.com/jmyu/CERIS-JGRA-Workshop  (commit 730af93), folder 1Sorghum:
#   7Envs_envParas_DAP122.txt  daily environmental parameters, 7 environments,
#                              122 days after planting
#   Trait_records.txt          flowering time of 237 recombinant inbred lines
# Clone or download that repository and pass the 1Sorghum folder as the first
# argument (default: ./CERIS-JGRA-Workshop/1Sorghum).
#
# Search, as in the workshop code (Exhaustive_search in Sub_functions_Workshop.r):
#   trait          flowering time in growing degree days (FTgdd), averaged over
#                  lines within each environment (the environment means)
#   parameters     day length (DL), growing degree days (GDD), photothermal time
#                  (PTT), photothermal ratio (PTR)
#   windows        every start day d1 and end day d2 with d2 >= d1 + 6 within
#                  days 1 to 122 after planting (6,786 windows per parameter)
#   candidate      mean of the parameter over days d1..d2 in each environment
#   rho            Pearson correlation of the candidate with the environment means
# 4 parameters x 6,786 windows = 27,144 candidates.
#
# Output: ./sorghum_scan/sorghum_ceris_scan.csv (para, d1, d2, rho, abs_rho,
#   tau_max) and sorghum_ceris_scan_summary.txt. Runs in a few seconds.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path("CERIS-JGRA-Workshop", "1Sorghum")
out_dir  <- if (length(args) >= 2) args[2] else "sorghum_scan"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

env_file   <- file.path(data_dir, "7Envs_envParas_DAP122.txt")
trait_file <- file.path(data_dir, "Trait_records.txt")
stopifnot(file.exists(env_file), file.exists(trait_file))

env_paras <- read.table(env_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
traits    <- read.table(trait_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        na.strings = "NA")

trait  <- "FTgdd"
paras  <- c("PTR", "DL", "GDD", "PTT")
n_days <- 122
min_span <- 6

env_means <- tapply(traits[[trait]], traits$env_code, mean, na.rm = TRUE)
envs <- names(env_means)
y <- as.numeric(env_means)

# days x environments matrix of one parameter, rows in day order after planting
para_matrix <- function(p) {
  sapply(envs, function(e) {
    v <- env_paras[env_paras$env_code == e, p]
    stopifnot(length(v) >= n_days)
    v[seq_len(n_days)]
  })
}

windows <- do.call(rbind, lapply(1:(n_days - min_span), function(d1)
  cbind(d1 = d1, d2 = (d1 + min_span):n_days)))

res <- do.call(rbind, lapply(paras, function(p) {
  M  <- para_matrix(p)
  cs <- rbind(0, apply(M, 2, cumsum))
  W  <- (cs[windows[, "d2"] + 1, , drop = FALSE] - cs[windows[, "d1"], , drop = FALSE]) /
        (windows[, "d2"] - windows[, "d1"] + 1)
  rho <- apply(W, 1, function(w) cor(w, y))
  data.frame(para = p, d1 = windows[, "d1"], d2 = windows[, "d2"], rho = rho)
}))
res$abs_rho <- abs(res$rho)
res$tau_max <- sqrt(pmax(0, 1 - res$rho^2))

write.csv(res, file.path(out_dir, "sorghum_ceris_scan.csv"), row.names = FALSE)

best <- res[which.max(res$abs_rho), ]
summ <- c(
  sprintf("Candidates: %d (%d parameters x %d windows)", nrow(res), length(paras), nrow(windows)),
  sprintf("Environments: %d; lines with records: %d", length(envs), length(unique(traits$line_code))),
  sprintf("Highest |rho|: %s, days %d to %d, rho = %.4f, tau_max = %.3f, 1/tau_max^2 = %.0f",
          best$para, best$d1, best$d2, best$rho, best$tau_max, 1 / best$tau_max^2),
  sprintf("Candidates with |rho| >= 0.99: %d; >= 0.90: %d",
          sum(res$abs_rho >= 0.99), sum(res$abs_rho >= 0.90)),
  paste("Candidates with |rho| >= 0.99 by parameter:",
        paste(names(table(res$para[res$abs_rho >= 0.99])), table(res$para[res$abs_rho >= 0.99]),
              sep = " ", collapse = "; "))
)
writeLines(summ, file.path(out_dir, "sorghum_ceris_scan_summary.txt"))
cat(summ, sep = "\n")
