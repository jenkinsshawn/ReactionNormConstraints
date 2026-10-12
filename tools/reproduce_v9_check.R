# =============================================================================
# tools/reproduce_v9_check.R  (ReactionNormConstraints v10.0.0)
#
# Checks that the v10 code reproduces v9.0.0 results with v10 options off.
#
#   Part 1: re-runs selected conditions of Experiment 1 (rho sweep) at full
#           replicate count (1000) with the v9.0.0 seeds and compares the
#           summary row with the committed sim_output_v9/summary_rho_sweep.csv.
#   Part 2 (optional): if RN_V9_DIR points to a folder of per-replicate v9.0.0
#           experiment tables (rho_sweep.csv and others), re-runs a sample of
#           those replicates (RN_REPRO_ROWS per table; 0 = all) with v10 and the
#           same seeds and compares every common column. Set SMOKE_TEST to match
#           the mode the v9 tables were produced in (it fixes the replicate
#           count and therefore the seed-to-condition map).
#
# Run from the repository root:  Rscript tools/reproduce_v9_check.R
# Env: RN_REPRO_RHO (comma list, default "0.996"; "none" skips Part 1), RN_V9_DIR,
#      RN_REPRO_ROWS (default 40),
#      RN_REPRO_OUT (default reproduce_check_output)
# =============================================================================

src <- readLines("sim_reaction_norm_gwas_v10.R", warn = FALSE)
kill <- grepl("^\\s*if\\s*\\(\\s*!\\s*interactive\\s*\\(\\s*\\)\\s*\\)\\s*main\\s*\\(\\s*\\)", src)
eval(parse(text = paste(src[!kill], collapse = "\n")), envir = globalenv())
OUT <- Sys.getenv("RN_REPRO_OUT", "reproduce_check_output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rho_check <- as.numeric(strsplit(Sys.getenv("RN_REPRO_RHO", "0.996"), ",")[[1]])
rho_check <- rho_check[!is.na(rho_check)]

# ---- Part 1 ----
if (length(rho_check)) {
committed <- read.csv(file.path("sim_output_v9", "summary_rho_sweep.csv"))
rows <- list(); t1 <- proc.time()
for (rv in rho_check) {
  i <- which(abs(rho_grid_fine - rv) < 1e-9)
  if (length(i) != 1L) stop("rho ", rv, " is not a single point of rho_grid_fine")
  n_full <- 1000L
  res <- run_experiment(data.frame(rho = rho_grid_fine[i]), default_params, n_full,
                        20260221L + (i - 1L) * n_full, sprintf("Exp 1 condition %d (rho = %s)", i, rv))
  new <- summarize_rho_sweep(res)
  old <- committed[abs(committed$rho_target - round(rv, 4)) < 1e-9, ]
  cols <- intersect(names(new), names(old))
  num  <- cols[vapply(cols, function(cl) is.numeric(new[[cl]]), logical(1))]
  d <- data.frame(rho = rv, column = num, committed_v9 = unlist(old[1, num]), v10 = unlist(new[1, num]),
                  stringsAsFactors = FALSE)
  d$abs_diff <- abs(d$v10 - d$committed_v9)
  rows[[length(rows) + 1]] <- d
}
p1 <- do.call(rbind, rows); rownames(p1) <- NULL
write.csv(p1, file.path(OUT, "part1_summary_row_comparison.csv"), row.names = FALSE)
cat(sprintf("Part 1: %d columns compared; max |diff| = %.3g; %.1f s\n", nrow(p1),
            max(p1$abs_diff, na.rm = TRUE), (proc.time() - t1)[3]))
}

# ---- Part 2 ----
# Seed-matched replicate identity. Every v9 replicate row carries its seed; the
# condition index follows from seed = offset + (i - 1) * n_reps + r. Each sampled
# v9 row is re-run with v10 run_replicate() under the same parameters and seed,
# and every common column is compared after the same CSV round trip.
v9dir <- Sys.getenv("RN_V9_DIR", Sys.getenv("RN_V9_SMOKE_DIR", ""))
if (nzchar(v9dir)) {
  n_rows <- suppressWarnings(as.integer(Sys.getenv("RN_REPRO_ROWS", "40")))
  pNL <- default_params; pNL$var_a_quad <- 0.03
  specs <- list(
    list(file = "rho_sweep.csv", design = data.frame(rho = rho_grid_fine), offset = 20260221L, n = rep_mult(1000), base = default_params),
    list(file = "kappa_comparison.csv", design = expand.grid(rho = rho_grid_coarse, kappa = c(0, 0.25, 0.5, 0.75, 1)), offset = 30000000L, n = rep_mult(200), base = default_params),
    list(file = "ds_coupling_sweep.csv", design = expand.grid(rho = rho_grid_coarse, ds_coupling = c(0, 0.3, 0.6, 0.9, 1)), offset = 50000000L, n = rep_mult(200), base = default_params),
    list(file = "pop_structure.csv", design = expand.grid(rho = rho_grid_coarse, pop_structure = c("diverse", "biparental", "nam")), offset = 65000000L, n = rep_mult(150), base = default_params),
    list(file = "nonlinear_test.csv", design = expand.grid(rho = rho_grid_coarse, nonlinear = c(FALSE, TRUE)), offset = 80000000L, n = rep_mult(200), base = pNL),
    list(file = "heterogeneous_residuals.csv", design = expand.grid(rho = rho_grid_sensitivity, hetero_res = c(FALSE, TRUE)), offset = 110000000L, n = rep_mult(100), base = default_params),
    list(file = "pca_lmm_gwas.csv", design = expand.grid(rho = rho_grid_sensitivity, gwas_method = c("simple", "pca", "lmm"), pop_structure = c("diverse", "nam")), offset = 115000000L, n = rep_mult(100), base = default_params),
    list(file = "missing_data.csv", design = expand.grid(rho = rho_grid_sensitivity, pct_missing = c(0, 0.05, 0.1, 0.2)), offset = 125000000L, n = rep_mult(100), base = default_params),
    list(file = "replication_sweep.csv", design = expand.grid(rho = rho_grid_sensitivity, n_reps = c(1L, 2L, 3L, 4L)), offset = 130000000L, n = rep_mult(100), base = default_params))
  set.seed(7L)
  jobs <- list(); meta <- list()
  for (sp in specs) {
    f <- file.path(v9dir, sp$file); if (!file.exists(f)) next
    old <- read.csv(f, stringsAsFactors = FALSE)
    take <- if (is.na(n_rows) || n_rows <= 0 || n_rows >= nrow(old)) seq_len(nrow(old)) else sort(sample(nrow(old), n_rows))
    for (r in take) {
      i <- (old$seed[r] - sp$offset - 1L) %/% sp$n + 1L
      p <- sp$base; for (nm in names(sp$design)) { v <- sp$design[[nm]][i]; p[[nm]] <- if (is.factor(v)) as.character(v) else v }
      jobs[[length(jobs) + 1]] <- list(params = p, seed = old$seed[r])
      meta[[length(meta) + 1]] <- list(file = sp$file, row = r)
    }
  }
  t2 <- proc.time()
  cl <- parallel::makeCluster(min(N_CORES, length(jobs)))
  parallel::clusterExport(cl, .worker_fns, envir = globalenv())
  new_rows <- parallel::parLapply(cl, jobs, function(j) run_replicate(j$params, j$seed))
  parallel::stopCluster(cl)
  out <- list()
  for (fl in unique(vapply(meta, `[[`, "", "file"))) {
    k <- which(vapply(meta, `[[`, "", "file") == fl)
    old <- read.csv(file.path(v9dir, fl), stringsAsFactors = FALSE)[vapply(meta[k], `[[`, 0L, "row"), ]
    new <- do.call(rbind, new_rows[k])
    new <- new[, !(names(new) %in% .v10_optional_cols), drop = FALSE]
    tmp <- tempfile(fileext = ".csv"); write.csv(new, tmp, row.names = FALSE)
    new <- read.csv(tmp, stringsAsFactors = FALSE)
    common <- intersect(names(old), names(new))
    same <- vapply(common, function(cl) isTRUE(all.equal(old[[cl]], new[[cl]], tolerance = 0, check.attributes = FALSE)), logical(1))
    num <- common[vapply(common, function(cl) is.numeric(old[[cl]]) && is.numeric(new[[cl]]), logical(1))]
    dif <- vapply(num, function(cl) { d <- abs(old[[cl]] - new[[cl]]); if (all(is.na(d))) 0 else max(d, na.rm = TRUE) }, numeric(1))
    nad <- vapply(common, function(cl) !identical(is.na(old[[cl]]), is.na(new[[cl]])), logical(1))
    disc <- c("n_sig_slope", "n_sig_int", "n_overlap", "n_unique_slope", grep("^power_|^fpr_|^fdr_power|^win_power", common, value = TRUE))
    disc <- intersect(disc, common)
    out[[fl]] <- data.frame(file = fl, rows_checked = length(k), rows_in_file = nrow(read.csv(file.path(v9dir, fl))),
                            cols_compared = length(common), cols_identical = sum(same),
                            detection_cols_identical = sprintf("%d/%d", sum(same[disc]), length(disc)),
                            max_abs_diff = max(dif), na_pattern_differs = sum(nad),
                            differing = paste(common[!same], collapse = " "), stringsAsFactors = FALSE)
  }
  p2 <- do.call(rbind, out); rownames(p2) <- NULL
  write.csv(p2, file.path(OUT, "part2_seed_matched_identity.csv"), row.names = FALSE)
  print(p2[, c("file", "rows_checked", "cols_compared", "cols_identical", "detection_cols_identical", "max_abs_diff")], row.names = FALSE)
  cat(sprintf("Part 2: %d replicates re-run in %.1f s\n", length(jobs), (proc.time() - t2)[3]))
}
