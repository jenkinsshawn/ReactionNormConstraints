# =============================================================================
# tools/derivation_check.R  (ReactionNormConstraints v10.0.0)
#
# Numerical check of the power derivation (Appendix A of the paper). Compares observed per-locus power
# from the simulation with (i) the conditional analytic prediction computed from
# the realized index and generating values in each replicate (sim Section 7C,
# analytic_power_conditional) and (ii) the closed-form expectation
# (closed_form_power). Neither prediction uses the GWAS outcome.
#
#   Part A: simplified core simulation (independent loci, MAF 0.30, fixed effect
#           magnitudes), grid over rho x mean-axis G x E x residual SD x
#           within-environment replicates x slope-only axis.
#   Part B: default (v9) simulation with LD, variable MAF and normal effects,
#           over rho, to check the conditional prediction in the realistic model.
#   Part C: per-locus variance explained by class at default parameters.
#
# Run from the repository root:  Rscript tools/derivation_check.R
# Env: RN_CHECK_REPS_A (default 30), RN_CHECK_REPS_B (default 40),
#      RN_CHECK_REPS_C (default 40), RN_CHECK_OUT (default derivation_check_output)
# =============================================================================

load_sim_quiet <- function(path = "sim_reaction_norm_gwas_v10.R") {
  src <- readLines(path, warn = FALSE)
  kill <- grepl("^\\s*if\\s*\\(\\s*!\\s*interactive\\s*\\(\\s*\\)\\s*\\)\\s*main\\s*\\(\\s*\\)", src)
  eval(parse(text = paste(src[!kill], collapse = "\n")), envir = globalenv())
}
load_sim_quiet()
geti <- function(v, d) { x <- suppressWarnings(as.integer(Sys.getenv(v, ""))); if (is.na(x)) d else x }
REPS_A <- geti("RN_CHECK_REPS_A", 30L); REPS_B <- geti("RN_CHECK_REPS_B", 40L); REPS_C <- geti("RN_CHECK_REPS_C", 40L)
OUT <- Sys.getenv("RN_CHECK_OUT", "derivation_check_output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t_all <- proc.time()

# ---------------- Part A ----------------
pA <- default_params; pA$core_sim <- TRUE; pA$fixed_effects <- TRUE
pA$n_so_qtl <- 5L; pA$predict_power <- TRUE; pA$run_gp <- FALSE
scen <- data.frame(gxe_mean_axis = c("default", "none"), ds_coupling = c(0.8, 0),
                   n_ap_qtl = c(2L, 0L), n_cn_qtl = c(3L, 0L), stringsAsFactors = FALSE)
g <- expand.grid(rho = c(0, 0.5, 0.85, 0.95, 0.99, 0.996), sigma_res = c(1, 3), n_reps = c(1L, 4L),
                 so_axis = c("index", "driver"), s = 1:2, stringsAsFactors = FALSE)
desA <- cbind(g[, c("rho", "sigma_res", "n_reps", "so_axis")], scen[g$s, ])
rownames(desA) <- NULL
tA <- proc.time()
resA <- run_experiment(desA, pA, REPS_A, 210000000L, "Derivation check A (core)")
tA <- (proc.time() - tA)[3]

cls_map <- c(indep = "power_indep_in_slope", shared = "power_shared_in_slope", so = "power_so_in_slope",
             cn = "power_cn_in_slope", ap = "power_ap_in_slope")
summ_cell <- function(d, p_cell, cf_random) {
  cf <- closed_form_power(p_cell, rho = d$rho_target[1], kappa = p_cell$kappa, random_effects = cf_random)
  rows <- list()
  for (cl in names(cls_map)) {
    if (!(cls_map[[cl]] %in% names(d)) || all(is.na(d[[cls_map[[cl]]]]))) next
    rows[[length(rows) + 1]] <- data.frame(test = "slope", cls = cl,
      obs = mean(d[[cls_map[[cl]]]], na.rm = TRUE),
      se = sd(d[[cls_map[[cl]]]], na.rm = TRUE) / sqrt(sum(!is.na(d[[cls_map[[cl]]]]))),
      pred_conditional = mean(d[[paste0("pred_power_", cl, "_slope")]], na.rm = TRUE),
      pred_closed_form = if (cl %in% cf$cls) cf$power_slope[cf$cls == cl] else NA_real_,
      ncp2_conditional = mean(d[[paste0("pred_ncp2_", cl, "_slope")]], na.rm = TRUE),
      ncp2_closed_form = if (cl %in% cf$cls) cf$ncp2_b[cf$cls == cl] else NA_real_,
      stringsAsFactors = FALSE)
  }
  for (cl in c("mean", "so")) {
    oc <- paste0(if (cl == "mean") "power_mean_in_int" else "power_so_in_int")
    if (!(oc %in% names(d))) next
    rows[[length(rows) + 1]] <- data.frame(test = "intercept", cls = cl,
      obs = mean(d[[oc]], na.rm = TRUE), se = sd(d[[oc]], na.rm = TRUE) / sqrt(nrow(d)),
      pred_conditional = mean(d[[paste0("pred_power_", cl, "_int")]], na.rm = TRUE),
      pred_closed_form = if (cl %in% cf$cls) cf$power_int[cf$cls == cl] else NA_real_,
      ncp2_conditional = NA_real_, ncp2_closed_form = if (cl %in% cf$cls) cf$ncp2_a[cf$cls == cl] else NA_real_,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$D_slope_pred_cond <- mean(d$pred_D_slope); out$D_slope_obs <- mean(d$obs_var_slope)
  out$D_slope_closed <- cf$D_slope[1]
  out
}
keyA <- interaction(resA$rho_target, resA$sigma_res, resA$n_reps, resA$so_axis, resA$gxe_mean_axis, drop = TRUE)
sumA <- do.call(rbind, lapply(split(resA, keyA), function(d) {
  p <- pA; for (nm in c("sigma_res", "n_reps", "so_axis", "ds_coupling", "n_ap_qtl", "n_cn_qtl")) p[[nm]] <- d[[nm]][1]
  cbind(rho = d$rho_target[1], tau = sqrt(1 - d$rho_target[1]^2), sigma_res = d$sigma_res[1], n_rep = d$n_reps[1],
        so_axis = d$so_axis[1], gxe_mean_axis = d$gxe_mean_axis[1], n_sim = nrow(d),
        summ_cell(d, p, cf_random = FALSE))
}))
rownames(sumA) <- NULL
write.csv(resA, file.path(OUT, "partA_core_replicates.csv"), row.names = FALSE)
write.csv(sumA, file.path(OUT, "partA_core_power_check.csv"), row.names = FALSE)

# ---------------- Part B ----------------
pB <- default_params; pB$predict_power <- TRUE; pB$run_gp <- FALSE
desB <- data.frame(rho = c(0, 0.5, 0.85, 0.95, 0.99, 0.996))
tB <- proc.time()
resB <- run_experiment(desB, pB, REPS_B, 220000000L, "Derivation check B (default sim)")
tB <- (proc.time() - tB)[3]
sumB <- do.call(rbind, lapply(split(resB, resB$rho_target), function(d)
  cbind(rho = d$rho_target[1], n_sim = nrow(d), summ_cell(d, pB, cf_random = TRUE))))
rownames(sumB) <- NULL
write.csv(sumB, file.path(OUT, "partB_default_power_check.csv"), row.names = FALSE)

# ---------------- Part C ----------------
pC <- default_params; pC$pve_cols <- TRUE; pC$run_gp <- FALSE
desC <- expand.grid(rho = c(0.5, 0.95, 0.996), n_so_qtl = c(0L, 5L))
tC <- proc.time()
resC <- run_experiment(desC, pC, REPS_C, 230000000L, "Derivation check C (variance explained)")
tC <- (proc.time() - tC)[3]
pve_cols <- grep("^pve_", names(resC), value = TRUE)
sumC <- do.call(rbind, lapply(split(resC, interaction(resC$n_so_qtl, resC$rho_target, drop = TRUE)), function(d) {
  r <- data.frame(n_so_qtl = d$n_so_qtl[1], rho = d$rho_target[1], n_sim = nrow(d))
  for (cl in pve_cols) r[[cl]] <- mean(d[[cl]], na.rm = TRUE)
  r }))
rownames(sumC) <- NULL
write.csv(sumC, file.path(OUT, "partC_locus_pve.csv"), row.names = FALSE)

# ---------------- accuracy summary ----------------
acc <- function(d, lab) {
  ok <- !is.na(d$pred_conditional)
  okc <- !is.na(d$pred_closed_form)
  data.frame(part = lab, n_cells = sum(ok),
             rmse_conditional = sqrt(mean((d$obs[ok] - d$pred_conditional[ok])^2)),
             max_abs_conditional = max(abs(d$obs[ok] - d$pred_conditional[ok])),
             r_conditional = cor(d$obs[ok], d$pred_conditional[ok]),
             rmse_closed_form = sqrt(mean((d$obs[okc] - d$pred_closed_form[okc])^2)),
             max_abs_closed_form = max(abs(d$obs[okc] - d$pred_closed_form[okc])),
             r_closed_form = cor(d$obs[okc], d$pred_closed_form[okc]))
}
accs <- rbind(acc(sumA, "A_core"), acc(sumB, "B_default"))
accs$seconds <- c(tA, tB)
write.csv(accs, file.path(OUT, "accuracy_summary.csv"), row.names = FALSE)
write.csv(data.frame(part = c("A", "B", "C", "total"), seconds = c(tA, tB, tC, (proc.time() - t_all)[3]),
                     reps = c(REPS_A * nrow(desA), REPS_B * nrow(desB), REPS_C * nrow(desC), NA),
                     workers = N_CORES), file.path(OUT, "timing.csv"), row.names = FALSE)
print(accs)
cat(sprintf("Total %.1f s\n", (proc.time() - t_all)[3]))
