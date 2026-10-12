# Closed-form power for independent-axis loci at the 126 design points of the kappa, panel-size and
# environment-number experiments, using default_params and closed_form_power() from the released
# simulation code (ReactionNormConstraints v10.0.0). Runtime about 45 s. Run from the repository root:
#   Rscript tools/closed_form_design_experiments.R
# Writes sim_output_v10/closed_form_design.csv with observed power beside the prediction.
ex <- parse("sim_reaction_norm_gwas_v10.R")
nm <- vapply(ex, function(e) if (is.call(e) && identical(e[[1]], as.name("<-"))) paste(deparse(e[[2]]), collapse = "") else "", "")
for (e in ex[nm %in% c("default_params", "closed_form_power")]) eval(e, globalenv())
P <- default_params
grid <- rbind(
  data.frame(exp = "kappa", expand.grid(rho = c(0, .5, .7, .85, .95, .99, .996), kappa = c(0, .25, .5, .75, 1)), n_geno = 500, n_env = 12),
  data.frame(exp = "panel", rho = rep(c(0, .5, .85, .95, .996), 7), kappa = 1,
             n_geno = rep(c(50, 100, 200, 500, 1000, 2000, 5000), each = 5), n_env = 12),
  data.frame(exp = "nenv", rho = rep(c(0, .5, .7, .85, .95, .99, .996), 8), kappa = 1, n_geno = 500,
             n_env = rep(c(3, 4, 5, 6, 8, 12, 16, 20), each = 7)))
grid$pred_closed_form <- sapply(seq_len(nrow(grid)), function(i) {
  p <- P; p$n_geno <- grid$n_geno[i]; p$n_env <- grid$n_env[i]
  cf <- suppressWarnings(closed_form_power(p, rho = grid$rho[i], kappa = grid$kappa[i]))
  cf$power_slope[cf$cls == "indep"]
})
obs <- function(f, by) {
  d <- read.csv(file.path("sim_output_v10", f))
  a <- aggregate(d$power_indep_in_slope, d[by], mean); names(a)[ncol(a)] <- "observed"; a
}
k <- obs("kappa_comparison.csv", c("rho_target", "kappa"))
g <- obs("genotype_sweep.csv", c("rho_target", "n_geno"))
v <- obs("nenv_sweep.csv", c("rho_target", "n_env"))
key <- function(x) round(x, 4)
grid$observed <- NA_real_
for (i in seq_len(nrow(grid))) {
  r <- key(grid$rho[i])
  hit <- switch(grid$exp[i],
    kappa = k$observed[key(k$rho_target) == r & key(k$kappa) == key(grid$kappa[i])],
    panel = g$observed[key(g$rho_target) == r & g$n_geno == grid$n_geno[i]],
    nenv  = v$observed[key(v$rho_target) == r & v$n_env == grid$n_env[i]])
  if (length(hit) == 1) grid$observed[i] <- hit
}
ok <- !is.na(grid$observed)
cat(sprintf("Cells matched: %d of %d; RMSE %.4f; max abs error %.4f; r %.4f\n", sum(ok), nrow(grid),
    sqrt(mean((grid$observed[ok] - grid$pred_closed_form[ok])^2)),
    max(abs(grid$observed[ok] - grid$pred_closed_form[ok])),
    cor(grid$observed[ok], grid$pred_closed_form[ok])))
write.csv(grid, file.path("sim_output_v10", "closed_form_design.csv"), row.names = FALSE)
