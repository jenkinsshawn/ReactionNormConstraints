# =============================================================================
# tools/benchmark_runtime.R  (ReactionNormConstraints v10.0.0)
#
# Measures per-replicate cost under full parallel load for the configurations
# the experiments use, then predicts full-scale wall time for every experiment
# from run_experiment()'s chunking rule (tasks of 10 replicates; at least 20
# chunks; one parLapply round per W tasks within a chunk).
#
# Run from the repository root:  Rscript tools/benchmark_runtime.R
# Env: RN_N_CORES (workers), RN_BENCH_OUT (default benchmark_output)
# Output: benchmark_output/per_replicate_cost.csv, runtime_estimate.csv
# =============================================================================

src <- readLines("sim_reaction_norm_gwas_v10.R", warn = FALSE)
kill <- grepl("^\\s*if\\s*\\(\\s*!\\s*interactive\\s*\\(\\s*\\)\\s*\\)\\s*main\\s*\\(\\s*\\)", src)
eval(parse(text = paste(src[!kill], collapse = "\n")), envir = globalenv())
OUT <- Sys.getenv("RN_BENCH_OUT", "benchmark_output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
W <- N_CORES
P <- function(...) { p <- default_params; a <- list(...); for (nm in names(a)) p[[nm]] <- a[[nm]]; p }
cfg <- list(
  default_gp   = P(), default_nogp = P(run_gp = FALSE),
  n1000 = P(n_geno = 1000), n2000 = P(n_geno = 2000), n5000 = P(n_geno = 5000),
  m5000 = P(n_markers = 5000), m10000 = P(n_markers = 10000),
  env20 = P(n_env = 20), biparental = P(pop_structure = "biparental"), nam = P(pop_structure = "nam"),
  lmm = P(gwas_method = "lmm"), pca = P(gwas_method = "pca"), missing20 = P(pct_missing = 0.2),
  reps4 = P(n_reps = 4L), qtl1000 = P(n_mean_qtl = 330, n_shared_qtl = 330, n_indep_qtl = 170, n_cn_qtl = 100, n_ap_qtl = 70),
  so5_pred = P(n_so_qtl = 5L, predict_power = TRUE),
  core_nogp = P(core_sim = TRUE, fixed_effects = TRUE, n_so_qtl = 5L, predict_power = TRUE, run_gp = FALSE),
  pve_nogp = P(pve_cols = TRUE, run_gp = FALSE), audit_nogp = P(audit_cols = TRUE, run_gp = FALSE))
cl <- parallel::makeCluster(W); on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c(.worker_fns, "rotation_test", "multi_index_gwas"), envir = globalenv())
rounds <- 2L; per_task <- 2L
cost <- data.frame(config = names(cfg), sec_per_rep = NA_real_)
for (k in seq_along(cfg)) {
  tasks <- lapply(seq_len(W * rounds), function(i) list(p = cfg[[k]], seeds = 5000L + (i - 1L) * per_task + seq_len(per_task)))
  t0 <- proc.time()[3]
  invisible(parallel::parLapply(cl, tasks, function(tk) for (s in tk$seeds) run_replicate(tk$p, s)))
  cost$sec_per_rep[k] <- (proc.time()[3] - t0) / (rounds * per_task)
  cat(sprintf("  %-13s %.3f s per replicate per worker under load\n", names(cfg)[k], cost$sec_per_rep[k]))
}
t0 <- proc.time()[3]
invisible(parallel::parLapply(cl, seq_len(W), function(i) {
  set.seed(i); p <- default_params; p$rho <- 0.95
  g <- simulate_genotypes(p$n_geno, p$n_markers, p$maf_min, p$n_mean_qtl, p$n_shared_qtl, p$n_indep_qtl, p$n_cn_qtl, p$n_ap_qtl, p$pop_structure)
  ph <- generate_phenotypes(g, p); ix <- construct_index(ph, 0.95, 1); rn <- estimate_reaction_norms(ph$Y, ix$x)
  rotation_test(g, rn, ix$actual_rho, 200L, p$alpha_bonf); multi_index_gwas(g, ph, ix, p$alpha_bonf) }))
cost <- rbind(cost, data.frame(config = "rotation_multi", sec_per_rep = proc.time()[3] - t0))
# single-worker (unloaded) cost of the default configuration, to model chunks
# that keep only a few workers busy
t0 <- proc.time()[3]
invisible(parallel::parLapply(cl, list(list(p = cfg$default_gp, seeds = 9001:9003)), function(tk) for (s in tk$seeds) run_replicate(tk$p, s)))
cost <- rbind(cost, data.frame(config = "default_gp_single_worker", sec_per_rep = (proc.time()[3] - t0) / 3))
t0 <- proc.time()[3]; cl2 <- parallel::makeCluster(W); parallel::clusterExport(cl2, .worker_fns, envir = globalenv()); parallel::stopCluster(cl2)
startup <- proc.time()[3] - t0
cost <- rbind(cost, data.frame(config = "cluster_startup", sec_per_rep = startup))
write.csv(cost, file.path(OUT, "per_replicate_cost.csv"), row.names = FALSE)
cc <- setNames(cost$sec_per_rep, cost$config)

# ---- predicted wall time per experiment (full scale) ----
# Per-replicate cost rises with the number of busy workers (shared cores and
# memory bandwidth). Costs above are at full load; a chunk that keeps only k of W
# workers busy runs at cost * (1/L + (1 - 1/L) * k / W), where L is the measured
# full-load / single-worker ratio for the default configuration.
L <- cc[["default_gp"]] / cc[["default_gp_single_worker"]]
exp_time <- function(n_cond, reps, t_rep, calls = 1L) {
  n_tasks <- n_cond * ceiling(reps / 10)
  n_chunks <- min(max(20L, ceiling(n_tasks / 50)), n_tasks)
  per_chunk <- ceiling(n_tasks / n_chunks)
  busy <- min(per_chunk, W)
  t_eff <- t_rep * (1 / L + (1 - 1 / L) * busy / W)
  calls * (n_chunks * ceiling(per_chunk / W) * min(reps, 10) * t_eff + cc[["cluster_startup"]])
}
nlev <- function(v) length(v)
lerp_n <- function(n) approx(c(500, 1000, 2000, 5000), c(cc[["default_gp"]], cc[["n1000"]], cc[["n2000"]], cc[["n5000"]]), xout = n, rule = 2)$y
lerp_m <- function(m) approx(c(2000, 5000, 10000), c(cc[["default_gp"]], cc[["m5000"]], cc[["m10000"]]), xout = m, rule = 2)$y
lerp_e <- function(e) approx(c(12, 20), c(cc[["default_gp"]], cc[["env20"]]), xout = e, rule = 2)$y
nf <- length(rho_grid_fine); nc <- length(rho_grid_coarse); ns <- length(rho_grid_sensitivity)
est <- rbind(
  data.frame(exp = "01 rho_sweep", sec = exp_time(nf, 1000, cc[["default_gp"]])),
  data.frame(exp = "02 null_model", sec = exp_time(nf, 300, cc[["default_gp"]])),
  data.frame(exp = "03 kappa", sec = exp_time(nc * 5, 200, cc[["default_gp"]])),
  data.frame(exp = "04 ds_coupling", sec = exp_time(nc * 5, 200, cc[["default_gp"]])),
  data.frame(exp = "05 env_range", sec = exp_time(nc * 2, 200, cc[["default_gp"]])),
  data.frame(exp = "06 nenv", sec = sum(sapply(c(3, 4, 5, 6, 8, 12, 16, 20), function(e) exp_time(nc, 200, lerp_e(e)))) - 7 * cc[["cluster_startup"]]),
  data.frame(exp = "07 pop_structure", sec = exp_time(nc * 3, 150, mean(c(cc[["default_gp"]], cc[["biparental"]], cc[["nam"]])))),
  data.frame(exp = "08 architecture", sec = exp_time(nc, 200, cc[["default_gp"]], calls = 5L)),
  data.frame(exp = "09 genotype", sec = sum(sapply(c(50, 100, 200, 500, 1000, 2000, 5000), function(n) exp_time(ns, 100, lerp_n(n)))) - 6 * cc[["cluster_startup"]]),
  data.frame(exp = "10 nonlinear", sec = exp_time(nc * 2, 200, cc[["default_gp"]])),
  data.frame(exp = "11 effect_size", sec = exp_time(5 * 4, 200, cc[["default_gp"]])),
  data.frame(exp = "12-13 rotation+multi_index", sec = (nc * 50 + nc * 150) / W * cc[["rotation_multi"]] / 2 + 2 * cc[["cluster_startup"]]),
  data.frame(exp = "14 marker", sec = sum(sapply(c(500, 1000, 2000, 5000, 10000), function(m) exp_time(ns, 100, lerp_m(m)))) - 4 * cc[["cluster_startup"]]),
  data.frame(exp = "15 heritability", sec = exp_time(ns * 8, 100, cc[["default_gp"]])),
  data.frame(exp = "16 qtl_count", sec = exp_time(ns, 100, mean(c(rep(cc[["default_gp"]], 5), cc[["qtl1000"]])), calls = 6L)),
  data.frame(exp = "17 hetero_res", sec = exp_time(ns * 2, 100, cc[["default_gp"]])),
  data.frame(exp = "18 pca_lmm", sec = exp_time(ns * 6, 100, mean(c(cc[["default_gp"]], cc[["pca"]], cc[["lmm"]], cc[["nam"]])))),
  data.frame(exp = "19 effect_dist", sec = exp_time(ns, 100, cc[["default_gp"]], calls = 2L)),
  data.frame(exp = "20 missing", sec = exp_time(ns * 4, 100, mean(c(cc[["default_gp"]], cc[["missing20"]])))),
  data.frame(exp = "21 replication", sec = exp_time(ns * 4, 100, mean(c(cc[["default_gp"]], cc[["reps4"]])))),
  data.frame(exp = "22 convergence", sec = exp_time(1, 300, cc[["default_gp"]], calls = 3L)),
  data.frame(exp = "23 seed", sec = exp_time(1, 100, cc[["default_gp"]], calls = 15L)),
  data.frame(exp = "24 ld_window", sec = exp_time(ns, 100, cc[["default_gp"]], calls = 4L)),
  data.frame(exp = "25 maf", sec = exp_time(ns * 5, 100, cc[["default_gp"]])),
  data.frame(exp = "26 slope_only_qtl", sec = exp_time(nc * 6, 200, cc[["so5_pred"]])),
  data.frame(exp = "27 core_simulation", sec = exp_time(nc, 200, cc[["core_nogp"]] + (cc[["default_gp"]] - cc[["default_nogp"]]))),
  data.frame(exp = "28 within_env_variance", sec = exp_time(ns * 16 + ns * 12, 100, cc[["core_nogp"]])),
  data.frame(exp = "29 locus_pve", sec = exp_time(ns * 2, 200, cc[["pve_nogp"]])),
  data.frame(exp = "audit", sec = exp_time(3, 300, cc[["audit_nogp"]]) + exp_time(4, 300, cc[["audit_nogp"]]) +
               exp_time(3, 300, cc[["audit_nogp"]], calls = 5L) + exp_time(nc * 10, 200, cc[["audit_nogp"]])),
  data.frame(exp = "check (A 200, B 300, C 200 reps)", sec = exp_time(96, 200, cc[["core_nogp"]]) +
               exp_time(6, 300, cc[["default_nogp"]] * 1.1) + exp_time(6, 200, cc[["pve_nogp"]])))
est$minutes <- round(est$sec / 60, 1)
est$part <- ifelse(grepl("^(0[1-9]|1[0-9]|2[0-5])", est$exp), "sim_v9",
            ifelse(grepl("^2[6-9]", est$exp), "sim_v10", ifelse(est$exp == "audit", "audit", "check")))
est$workers <- W
write.csv(est, file.path(OUT, "runtime_estimate.csv"), row.names = FALSE)
print(est[, c("exp", "minutes")], row.names = FALSE)
print(aggregate(minutes ~ part, est, sum))
