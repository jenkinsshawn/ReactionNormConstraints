# =============================================================================
# reaction_norm_audit.R  (ReactionNormConstraints v10.0.0)
#
# Loads sim_reaction_norm_gwas_v10.R from the repository root, runs the audit
# stages, writes CSVs and one figure to ./audit_output_v10/, prints the tables.
#
# Run from the repository root:
#   Rscript reaction_norm_audit.R
#   SMOKE_TEST=TRUE Rscript reaction_norm_audit.R      (replicates scaled to ~2%)
# PowerShell:  $env:SMOKE_TEST = "TRUE"; Rscript reaction_norm_audit.R
#
# The simulation file is looked up in the working directory and next to this
# script; output goes to ./audit_output_v10/. The simulation's pinned RNGkind and
# SMOKE_TEST setting apply here, and the replicate counts below pass through its
# rep_mult(). The loader drops the simulation's final if (!interactive()) main()
# line before evaluating it, so loading does not start the full simulation; the
# file on disk is not modified. The extra columns used here (H2_slope_bx and
# related) come from run_replicate() with params$audit_cols = TRUE.
#
# Stages:
#   1 diagnostics  cor(slope, intercept) by population structure at rho = 0.996;
#                  intercept against genotype mean; numbers of significant
#                  slope and intercept markers; slope heritability
#   2 he_test      five QTL architectures at rho = 0.95, 0.99, 0.996:
#                  cor(slope, mean), overlap and independent-axis power
#   3 ds_sweep     differential-sensitivity coupling grid including the
#                  default 0.8
#   4 figure       genome-wide -log10(P) correlation against rho, simulation
#                  and published values (Tibbs-Cortes et al. 2024)
# =============================================================================


# =============================================================================
# CONFIG
# =============================================================================

SIM_FILE <- "sim_reaction_norm_gwas_v10.R"
OUTDIR   <- "audit_output_v10"              # relative to the repository root

options(error = NULL)                       # no Browse[1]> on stop()

STAGES <- c("diagnostics", "he_test", "ds_sweep", "figure")
.st <- Sys.getenv("RN_AUDIT_STAGES", "")
if (nzchar(.st)) STAGES <- trimws(strsplit(.st, ",")[[1]])

# Full-run replicate counts. They are passed through rep_mult() after the sim is
# loaded, so SMOKE_TEST=TRUE scales them to ~2% (minimum 2).
N_REPS_DIAG_FULL <- 300    # stage 1:  3x300 + 4x300      =  2,100 runs
N_REPS_HE_FULL   <- 300    # stage 2:  5 arch x 3 rho x300 =  4,500 runs
N_REPS_DS_FULL   <- 200    # stage 3:  7 rho x 10 cds x200 = 14,000 runs

RUN_GP <- FALSE       # genomic prediction CV adds most of the runtime and is not
# used by these stages. Set TRUE to add the gp_* columns.


# =============================================================================
# LOADER
# =============================================================================

script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = FALSE)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/")))
  getwd()
}

find_sim <- function(path = SIM_FILE) {
  cand <- unique(c(path, file.path(getwd(), basename(path)), file.path(script_dir(), basename(path))))
  hit  <- cand[file.exists(cand) & !dir.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
  stop("Cannot find '", basename(path), "'. Run this script from the repository root ",
       "or set SIM_FILE at the top of this script.")
}

load_sim <- function(path = SIM_FILE) {
  SIM_PATH <<- find_sim(path)
  cat("  Loader: reading", SIM_PATH, "\n")
  
  if (!grepl("^([A-Za-z]:|/|\\\\)", OUTDIR)) OUTDIR <<- file.path(dirname(SIM_PATH), OUTDIR)
  
  src <- readLines(SIM_PATH, warn = FALSE)
  src <- gsub("\r$", "", src)
  
  kill <- grepl("^\\s*setwd\\s*\\(", src) |
    grepl("^\\s*if\\s*\\(\\s*!\\s*interactive\\s*\\(\\s*\\)\\s*\\)\\s*main\\s*\\(\\s*\\)", src)
  if (any(kill)) {
    cat(sprintf("  Loader: stripped %d line(s) before evaluating:\n", sum(kill)))
    for (i in which(kill)) cat(sprintf("    line %d: %s\n", i, trimws(src[i])))
  }
  
  eval(parse(text = paste(src[!kill], collapse = "\n")), envir = globalenv())
  
  need <- c("default_params", "run_experiment", "run_replicate", "simulate_genotypes",
            "generate_phenotypes", "construct_index", "estimate_reaction_norms",
            "fast_gwas", "compute_metrics", "compute_qg_metrics",
            "genomic_prediction_cv", "estimate_ammi",
            "N_CORES", "rho_grid_coarse", ".worker_fns",
            "rep_mult", "SMOKE", "v10_optional_metrics")
  miss <- need[!vapply(need, exists, logical(1), envir = globalenv())]
  if (length(miss)) {
    stop("Loaded '", SIM_PATH, "' but these are still undefined: ",
         paste(miss, collapse = ", "),
         "\n  The file is probably not the version this script was written against.")
  }
  cat(sprintf("  Loader: OK. RNGkind = %s; SMOKE_TEST = %s; %d workers. Output -> %s\n\n",
              paste(RNGkind(), collapse = "/"), SMOKE, N_CORES, OUTDIR))
  invisible(TRUE)
}


# =============================================================================
# HELPERS
# =============================================================================

# Column means by group, ignoring NA within each column. The formula interface
# of aggregate() would drop a whole row whenever any one column is NA
# (overlap_prop is NA when a replicate has no significant slope markers).
agg_mean <- function(d, cols, by) {
  out <- do.call(rbind, lapply(split(d, d[[by]], drop = TRUE), function(g) {
    row <- data.frame(by_val = as.character(g[[by]][1]), n_reps = nrow(g),
                      stringsAsFactors = FALSE)
    names(row)[1] <- by
    for (cl in cols) {
      row[[cl]]                 <- round(mean(g[[cl]], na.rm = TRUE), 4)
      row[[paste0(cl, "_nNA")]] <- sum(is.na(g[[cl]]))
    }
    row
  }))
  rownames(out) <- NULL
  out
}

boot_ci <- function(v, n_boot = 2000) {
  v <- v[!is.na(v)]
  if (length(v) < 2) return(c(NA_real_, NA_real_, NA_real_))
  b <- replicate(n_boot, mean(sample(v, length(v), replace = TRUE)))
  c(mean(v), unname(quantile(b, 0.025)), unname(quantile(b, 0.975)))
}

# Mean and bootstrap 95% interval of yvar at each value of xvar.
agg_ci <- function(d, xvar, yvar, n_boot = 2000) {
  out <- do.call(rbind, lapply(split(d, d[[xvar]], drop = TRUE), function(g) {
    ci <- boot_ci(g[[yvar]], n_boot)
    data.frame(x = g[[xvar]][1], mn = ci[1], lo = ci[2], hi = ci[3])
  }))
  out <- out[order(out$x), ]
  rownames(out) <- NULL
  out
}

save_csv <- function(d, name) {
  if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
  f <- file.path(OUTDIR, name)
  write.csv(d, f, row.names = FALSE)
  cat(sprintf("  -> %s\n", f))
  invisible(f)
}

banner <- function(txt) {
  cat("\n", strrep("=", 78), "\n  ", txt, "\n", strrep("=", 78), "\n\n", sep = "")
}


# =============================================================================
# Stage 1 columns
#
# H2_slope is cor(true_slope_h, estimated slope)^2, where true_slope_h is the
# genotype's slope on the driver h. The index is x = rho*em_std + tau*z_orth, so
# the slope on x estimates approximately rho*true_slope_h + tau*true_slope_z.
# H2_slope_bx scores the estimate against that combination, and H2_slope_x
# against the exact genotypic slope on x, sum_j c_j G_ij with
# c_j = (x_j - mean(x)) / S_xx. n_sig_ratio_int_slope is the ratio of significant
# intercept to slope markers. These columns come from run_replicate() when
# params$audit_cols = TRUE (simulation Section 7C).
# =============================================================================


# =============================================================================
# Differential-sensitivity coupling grid
#
# The main experiment uses c_DS = 0, 0.3, 0.6, 0.9, 1. This grid adds 0.1 to 0.5
# and the default value 0.8.
# =============================================================================

run_ds_coupling_sweep_audit <- function(n_reps = N_REPS_DS, base = NULL) {
  p <- if (is.null(base)) default_params else base
  run_experiment(
    expand.grid(rho = rho_grid_coarse,
                ds_coupling = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0)),
    p, n_reps, 50000000L, "DS coupling sweep")
}

# overlap_prop is n_overlap / n_sig_slope. As c_DS rises, slope power rises and
# the added slope hits are small-effect differential-sensitivity loci that the
# intercept test also misses, so overlap_prop can fall while power rises.
# n_sig_slope is reported beside overlap_prop for that reason.
summarize_ds_sweep <- function(res, rho_at = 0.996) {
  d <- res[abs(res$rho_target - rho_at) < 1e-6, ]
  out <- do.call(rbind, lapply(split(d, d$ds_coupling), function(g) {
    data.frame(
      ds_coupling        = g$ds_coupling[1],
      n_reps             = nrow(g),
      overlap_prop       = round(mean(g$overlap_prop, na.rm = TRUE), 4),
      n_sig_slope        = round(mean(g$n_sig_slope, na.rm = TRUE), 2),
      n_sig_int          = round(mean(g$n_sig_int, na.rm = TRUE), 2),
      power_shared_slope = round(mean(g$power_shared_in_slope, na.rm = TRUE), 4),
      power_indep_slope  = round(mean(g$power_indep_in_slope, na.rm = TRUE), 4),
      cor_slope_int      = round(mean(g$cor_slope_int, na.rm = TRUE), 4)
    )
  }))
  out <- out[order(out$ds_coupling), ]; rownames(out) <- NULL
  out
}


# =============================================================================
# Five QTL architectures
#
# He et al. (2025) report phenotypic cor(slope, mean) of 0.00 to 0.20 at
# |rho| = 0.95 to 0.99, with slope-intercept marker Jaccard of 0.28 to 0.39.
# Effect sizes here are drawn as rnorm(nq, 0, sqrt(vq / (2*f*(1-f)))), which fixes
# each locus's variance contribution regardless of allele frequency, so changing
# the population design changes cor_slope_int little. The architectures below
# vary the numbers of mean-only, differential-sensitivity and antagonistic-
# pleiotropy loci and the coupling, and report cor(slope, mean), overlap and
# independent-axis power at each rho.
# =============================================================================

run_architecture_he_test <- function(n_reps = N_REPS_HE,
                                     rho_levels = c(0.95, 0.99, 0.996),
                                     base = NULL) {
  base_p <- if (is.null(base)) default_params else base
  configs <- list(
    default       = list(n_mean_qtl = 10, n_shared_qtl = 10, n_ap_qtl = 2,  ds_coupling = 0.8),
    ap_rich       = list(n_mean_qtl = 10, n_shared_qtl = 10, n_ap_qtl = 10, ds_coupling = 0.8),
    ds_poor       = list(n_mean_qtl = 10, n_shared_qtl = 3,  n_ap_qtl = 2,  ds_coupling = 0.8),
    mean_dominant = list(n_mean_qtl = 40, n_shared_qtl = 10, n_ap_qtl = 2,  ds_coupling = 0.8),
    weak_coupling = list(n_mean_qtl = 10, n_shared_qtl = 10, n_ap_qtl = 2,  ds_coupling = 0.2)
  )
  out <- vector("list", length(configs))
  for (ci in seq_along(configs)) {
    p <- base_p
    for (nm in names(configs[[ci]])) p[[nm]] <- configs[[ci]][[nm]]
    r <- run_experiment(data.frame(rho = rho_levels), p, n_reps,
                        140000000L + ci * 100000L,
                        paste("He-test:", names(configs)[ci]))
    r$architecture <- names(configs)[ci]
    out[[ci]] <- r
  }
  do.call(rbind, out)
}

summarize_he_test <- function(res) {
  key <- paste(res$architecture, res$rho_target, sep = "|")
  out <- do.call(rbind, lapply(split(res, key, drop = TRUE), function(d) {
    csm <- boot_ci(d$cor_slope_mean); ov <- boot_ci(d$overlap_prop)
    data.frame(
      architecture   = d$architecture[1],
      rho            = d$rho_target[1],
      cor_slope_mean = round(csm[1], 3), csm_lo = round(csm[2], 3), csm_hi = round(csm[3], 3),
      overlap_prop   = round(ov[1], 3),  ov_lo  = round(ov[2], 3),  ov_hi  = round(ov[3], 3),
      power_indep    = round(mean(d$power_indep_in_slope, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$rho, out$architecture), ]; rownames(out) <- NULL
  out
}


# =============================================================================
# Numbers of significant slope and intercept markers
# =============================================================================

summarize_power_asymmetry <- function(res) {
  out <- do.call(rbind, lapply(split(res, res$rho_target), function(g) {
    data.frame(
      rho                   = g$rho_target[1],
      n_reps                = nrow(g),
      mean_n_sig_int        = round(mean(g$n_sig_int, na.rm = TRUE), 2),
      mean_n_sig_slope      = round(mean(g$n_sig_slope, na.rm = TRUE), 2),
      ratio_int_slope       = round(mean(g$n_sig_int, na.rm = TRUE) /
                                      mean(g$n_sig_slope, na.rm = TRUE), 3),
      pct_reps_int_gt_slope = round(mean(g$n_sig_int > g$n_sig_slope, na.rm = TRUE), 3)
    )
  }))
  out <- out[order(out$rho), ]; rownames(out) <- NULL
  out
}


# =============================================================================
# Simulated curve with published values
#
# Plots the simulation mean and 95% interval of a metric against rho, with
# published values as points. Each point is the value the study reports; a dotted
# segment joins it to the simulated mean at the same rho.
# =============================================================================

fig_constraint_curve <- function(res_rho, emp = NULL, metric = "cor_slope_int",
                                 ylab = "Correlation (slope, intercept)",
                                 outfile = NULL) {
  if (!is.null(emp)) {
    if (!all(c("rho", "y", "lb", "tr") %in% names(emp)))
      stop("emp needs columns rho, y, lb, tr, where y is the value the study reports.")
    if (any(is.na(emp$y)))
      stop("emp$y is NA for: ", paste(emp$lb[is.na(emp$y)], collapse = ", "),
           ". Drop those rows or supply the published values.")
  }
  
  if (!is.null(outfile)) {
    if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
    tiff(outfile, width = 7.5, height = 5.5, units = "in", res = 300, compression = "lzw")
    on.exit(dev.off(), add = TRUE)
  }
  op <- par(mar = c(4.8, 5.2, 3.0, 1.5)); on.exit(par(op), add = TRUE)
  
  a1 <- agg_ci(res_rho, "rho_target", metric)
  yr <- range(c(a1$lo, a1$hi, if (!is.null(emp)) emp$y), na.rm = TRUE)
  yr <- c(min(0, yr[1]) - 0.02, max(yr[2] * 1.10, 0.2))
  
  plot(a1$x, a1$mn, type = "n", xlim = c(-0.03, 1.01), ylim = yr,
       xlab = expression("|" * rho * "|"), ylab = ylab,
       main = "Simulation and published values")
  rect(-0.05, yr[1], 0.85, yr[2], col = adjustcolor("#0072B2", 0.06), border = NA)
  rect( 0.85, yr[1], 0.95, yr[2], col = adjustcolor("#E69F00", 0.10), border = NA)
  rect( 0.95, yr[1], 1.05, yr[2], col = adjustcolor("#CC79A7", 0.10), border = NA)
  polygon(c(a1$x, rev(a1$x)), c(a1$lo, rev(a1$hi)),
          col = adjustcolor("#0072B2", 0.20), border = NA)
  lines(a1$x, a1$mn, lwd = 2.5, col = "#0072B2")
  
  if (!is.null(emp)) {
    pch_tr <- c(FT = 21, PH = 22, Multi = 23)
    col_tr <- c(FT = "#E69F00", PH = "#009E73", Multi = "#CC79A7")
    sim_at_emp <- approx(a1$x, a1$mn, xout = emp$rho, rule = 2)$y
    # Dotted segments join each published value to the simulated mean.
    segments(emp$rho, emp$y, emp$rho, sim_at_emp, lty = 3, col = adjustcolor("grey30", 0.8))
    points(emp$rho, emp$y, pch = pch_tr[emp$tr], bg = col_tr[emp$tr], col = "black", cex = 1.6)
    legend("topleft",
           legend = c("Simulation (95% CI)", "Observed: flowering", "Observed: height",
                      "Observed: multiple", "Observed - simulation"),
           pch = c(NA, 21, 22, 23, NA),
           pt.bg = c(NA, col_tr[["FT"]], col_tr[["PH"]], col_tr[["Multi"]], NA),
           col = c("#0072B2", "black", "black", "black", "grey30"),
           lwd = c(2.5, NA, NA, NA, 1), lty = c(1, NA, NA, NA, 3),
           bty = "n", pt.cex = 1.2, cex = 0.85)
    r_emp <- cor(sim_at_emp, emp$y, use = "complete.obs")
    mtext(sprintf("observed vs simulation: r = %.2f (n = %d)", r_emp, nrow(emp)),
          side = 3, line = 0.2, adj = 1, cex = 0.8, col = "grey30")
    cat(sprintf("\n  Figure: observed vs simulation r = %.2f over %d points.\n",
                r_emp, nrow(emp)))
  } else {
    legend("topleft", legend = "Simulation (95% CI)", col = "#0072B2",
           lwd = 2.5, lty = 1, bty = "n", cex = 0.85)
  }
  invisible(a1)
}

# Genome-wide -log10(P) correlation between slope and intercept GWAS reported by
# Tibbs-Cortes et al. (2024) for 13 maize nested association mapping traits.
tc24_logp_observed <- data.frame(
  rho = c(0.999, 0.998, 0.996, 0.993, 0.988, 0.998, 0.998, 0.993, 0.960, 0.980, 0.975, 0.999, 0.998),
  y   = c(0.11,  0.10,  0.04,  0.09,  0.09,  0.13,  0.14,  0.03,  0.02,  0.07,  0.21,  0.19,  0.08),
  lb  = c("DTA","DTS","ASI","EH","PH","GDD.DTA","GDD.DTS","GDD.ASI","LL","LW","TPBN","TKW","T20KW"),
  tr  = c("FT","FT","FT","PH","PH","FT","FT","FT","Multi","Multi","Multi","Multi","Multi"),
  stringsAsFactors = FALSE
)


# =============================================================================
# STAGES
# =============================================================================

stage_diagnostics <- function(p0, n_reps) {
  banner("STAGE 1: diagnostics")
  
  cat("-- 1a. cor(slope, intercept) at rho = 0.996 by population structure\n")
  cat("   (Biparental and NAM designs are not calibrated under the simple test.)\n\n")
  ps <- run_experiment(expand.grid(rho = 0.996,
                                   pop_structure = c("diverse", "biparental", "nam"),
                                   stringsAsFactors = FALSE),
                       p0, n_reps, 65000000L, "Pop structure at rho = 0.996")
  t1 <- agg_mean(ps, c("cor_slope_int", "cor_slope_mean", "overlap_prop"), "pop_structure")
  print(t1); save_csv(t1, "01_pop_structure.csv")
  
  cat("\n-- 1b. Estimated intercept against genotype mean\n")
  cat("   With kappa = 1, mean(x) = 0, so intercept = genotype mean.\n\n")
  g  <- simulate_genotypes(500, 2000, 0.05, 10, 10, 5, 3, 2, "diverse")
  ph <- generate_phenotypes(g, p0)
  ix <- construct_index(ph, 0.996, 1.0)
  rn <- estimate_reaction_norms(ph$Y, ix$x)
  cat(sprintf("   mean(x) = %.3e   max|intercept - mean_perf| = %.3e\n",
              mean(ix$x), max(abs(rn$intercept - rn$mean_perf))))
  cat(sprintf("   cor(slope, intercept) = %.6f   cor(slope, mean_perf) = %.6f\n",
              cor(rn$slope, rn$intercept), cor(rn$slope, rn$mean_perf)))
  
  cat("\n-- 1c. Numbers of significant intercept and slope markers\n\n")
  rs <- run_experiment(data.frame(rho = c(0.5, 0.85, 0.95, 0.996)),
                       p0, n_reps, 20260221L, "Power asymmetry")
  t3 <- summarize_power_asymmetry(rs)
  print(t3); save_csv(t3, "02_power_asymmetry.csv")
  
  cat("\n-- 1d. Slope heritability against two targets\n")
  cat("   H2_slope: against the slope on h. H2_slope_bx: against\n")
  cat("   rho*true_slope_h + tau*true_slope_z.\n\n")
  t4 <- agg_mean(rs, c("H2_int", "H2_slope", "H2_slope_bx", "H2_ratio", "H2_ratio_bx"),
                 "rho_target")
  print(t4); save_csv(t4, "03_heritability_targets.csv")
  
  save_csv(rs, "04_rho_sweep_raw.csv")
  invisible(list(pop_structure = ps, rho_sweep = rs))
}

stage_he_test <- function(p0, n_reps) {
  banner("STAGE 2: five QTL architectures")
  ht <- run_architecture_he_test(n_reps = n_reps, base = p0)
  s <- summarize_he_test(ht)
  print(s); save_csv(s, "05_he_test_summary.csv"); save_csv(ht, "06_he_test_raw.csv")
  invisible(ht)
}

stage_ds_sweep <- function(p0, n_reps) {
  banner("STAGE 3: differential-sensitivity coupling grid")
  ds <- run_ds_coupling_sweep_audit(n_reps = n_reps, base = p0)
  s <- summarize_ds_sweep(ds, rho_at = 0.996)
  print(s); save_csv(s, "07_ds_sweep_summary.csv"); save_csv(ds, "08_ds_sweep_raw.csv")
  invisible(ds)
}

stage_figure <- function(rs) {
  banner("STAGE 4: -log10(P) correlation, simulation and published values")
  f <- file.path(OUTDIR, "Figure_4_constraint_curve.tiff")
  fig_constraint_curve(rs, emp = tc24_logp_observed, metric = "cor_logp",
                       ylab = expression("Genome-wide " * italic(r) * " (-log"[10] * italic(P) * ")"),
                       outfile = f)
  cat(sprintf("  -> %s\n", f))
  cat("  Published values: Tibbs-Cortes et al. (2024), 13 maize traits.\n")
  invisible(f)
}


# =============================================================================
# MAIN
# =============================================================================

run_audit <- function() {
  t0 <- proc.time()
  banner("REACTION NORM AUDIT (v10.0.0)")
  
  load_sim()
  N_REPS_DIAG <<- rep_mult(N_REPS_DIAG_FULL)
  N_REPS_HE   <<- rep_mult(N_REPS_HE_FULL)
  N_REPS_DS   <<- rep_mult(N_REPS_DS_FULL)
  cat(sprintf("  Replicates per cell: diagnostics %d, he_test %d, ds_sweep %d\n",
              N_REPS_DIAG, N_REPS_HE, N_REPS_DS))
  
  p0 <- default_params
  p0$run_gp <- isTRUE(RUN_GP)
  p0$audit_cols <- TRUE
  
  # Check that the extra columns are produced before the full run
  chk <- run_replicate(p0, 1L)
  if (!all(c("H2_slope_bx", "H2_slope_x", "n_sig_ratio_int_slope") %in% names(chk)))
    stop("run_replicate did not return the audit columns.")
  cat("  Audit-column check: OK (H2_slope_bx, H2_slope_x present).\n")
  
  timing <- data.frame(stage = character(0), seconds = numeric(0), stringsAsFactors = FALSE)
  tick <- function(lbl, t) timing[nrow(timing) + 1L, ] <<- list(lbl, round((proc.time() - t)[3], 1))
  out <- list()
  if ("diagnostics" %in% STAGES) { t <- proc.time(); out$diag <- stage_diagnostics(p0, N_REPS_DIAG); tick("diagnostics", t) }
  if ("he_test"     %in% STAGES) { t <- proc.time(); out$he   <- stage_he_test(p0, N_REPS_HE);       tick("he_test", t) }
  if ("ds_sweep"    %in% STAGES) { t <- proc.time(); out$ds   <- stage_ds_sweep(p0, N_REPS_DS);      tick("ds_sweep", t) }
  if ("figure"      %in% STAGES && !is.null(out$diag)) {
    t <- proc.time(); out$fig <- stage_figure(out$diag$rho_sweep); tick("figure", t) }
  save_csv(timing, "00_audit_timing.csv")
  
  banner(sprintf("DONE in %.1f min. Output in %s",
                 (proc.time() - t0)[3] / 60, normalizePath(OUTDIR, mustWork = FALSE)))
  invisible(out)
}

if (!identical(Sys.getenv("RN_AUDIT_NO_RUN"), "TRUE")) audit <- run_audit()
