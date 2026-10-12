#!/usr/bin/env Rscript

# ============================================================================
# ReactionNormConstraints v10.0.0
# Simulation framework: structural constraints on environmental index-based
#   reaction norm GWAS
#
# Jenkins SD, Graef GL
# "Detection power for plasticity loci when the environmental index is chosen
#  to match environment means"
#
# v10.0.0 (additive release; see CHANGES.md):
#   - No-intercept slope QTL class (n_so_qtl > 0): loci whose effect changes the
#     slope along the index x (or along the environmental driver h) with zero
#     effect on the intercept at x = 0 or, separately, zero effect on the
#     environment-averaged genotype mean. Effects and positions are drawn from an
#     independent RNG substream, so switching the class on leaves every other
#     locus, effect, environment and residual of a replicate unchanged.
#   - Optional simplified core simulation (core_sim = TRUE): independent loci,
#     fixed allele frequency, fixed effect magnitudes (fixed_effects = TRUE), no LD.
#   - Optional parameter-free analytic power prediction per replicate
#     (predict_power = TRUE), from the slope and intercept contrasts of the
#     realized index and the true generating values (Appendix A of the paper).
#   - Optional per-locus variance explained by QTL class (pve_cols = TRUE).
#   - Optional audit columns used by reaction_norm_audit.R (audit_cols = TRUE).
#   - Experiments 26-29 (v10 additions) run after Experiments 1-25.
#   All additions are off at default_params; Experiments 1-25 reproduce v9.0.0
#   replicate-for-replicate (same seeds, same columns, same values).
#
# Dependencies: base R only. The 'parallel' package ships with R. Requires R >= 3.6.0.
# Output:       written to a relative ./sim_output_v10/ directory created at run time.
#               Run from the repository root; no setwd() is used.
#
# Reproducibility:
#   - The RNG regime is pinned in Section 0 so results reproduce exactly across R
#     installations (R 3.6.0 changed the default sample() algorithm).
#   - Per-replicate seeds are deterministic, so results do not depend on core count.
#   - sessionInfo() is written to sim_output_v10/run_metadata_and_session.txt.
#
# Run control (environment variables):
#   SMOKE_TEST=TRUE        scale all replicate/permutation/bootstrap counts to ~2%.
#   RN_RUN_V9_SET=FALSE    skip Experiments 1-25 (and their summaries/figures).
#   RN_RUN_V10_SET=FALSE   skip Experiments 26-29.
#   RN_N_CORES=<k>         override the number of parallel workers.
#
# Runtime: Experiments 1-25 take ~2-3 hours on a 12-core / 128 GB machine.
#   See run_full_v10.R for measured estimates of the full v10 run.
# ============================================================================

cat("
================================================================================
  Reaction Norm GWAS Simulation Framework v10.0.0
  Jenkins & Graef (2026)
  Structural Constraints on Environmental Index-Based Genetic Dissection
================================================================================
\n")

# ============================================================================
# SECTION 0: REPRODUCIBILITY SETUP
# ============================================================================

# Output is written to a relative ./sim_output_v10/ directory, created in main().
# Run this script from the repository root. No setwd() is used, so it runs
# unmodified on any machine.

# Pin the RNG regime. R 3.6.0 changed the default sample() algorithm; pinning it
# guarantees identical draws (and identical results) across R >= 3.6.0.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Smoke test: set SMOKE_TEST=TRUE to scale every replicate, permutation,
# and bootstrap count to ~2% for a fast end-to-end pipeline check. Leave unset for
# the full production sweep.
SMOKE <- tryCatch(as.logical(Sys.getenv("SMOKE_TEST", "FALSE")), error = function(e) FALSE)
if (length(SMOKE) != 1L || is.na(SMOKE)) SMOKE <- FALSE
rep_mult <- function(n) if (SMOKE) max(2L, as.integer(round(n * 0.02))) else as.integer(n)
if (SMOKE) cat("  [SMOKE_TEST] Reduced-replicate mode: counts scaled to ~2%.\n\n")

# ============================================================================
# SECTION 1: DEFAULT PARAMETERS & PARALLELIZATION
# ============================================================================

total_cores <- parallel::detectCores(logical = TRUE)
N_CORES <- max(1L, floor(total_cores * 0.85))
.env_cores <- suppressWarnings(as.integer(Sys.getenv("RN_N_CORES", "")))
if (length(.env_cores) == 1L && !is.na(.env_cores) && .env_cores >= 1L) N_CORES <- .env_cores

cat(sprintf("  Detected %d logical threads. Using %d workers for load balancing.\n\n", 
            total_cores, N_CORES))

default_params <- list(
  # --- Genotype dimensions ---
  n_geno       = 500,
  n_markers    = 2000,
  maf_min      = 0.05,
  pop_structure = "diverse",     
  
  # --- QTL architecture (five QTL classes, Des Marais et al. 2013) ---
  n_mean_qtl   = 10,             
  n_shared_qtl = 10,             
  n_indep_qtl  = 5,              
  n_cn_qtl     = 3,              
  n_ap_qtl     = 2,              
  
  # --- Effect size variances ---
  var_a_mean   = 0.04,
  var_a_shared = 0.04,
  ds_coupling  = 0.8,            
  var_a_indep  = 0.04,
  var_a_cn     = 0.04,
  cn_threshold = 0.0,            
  var_a_ap     = 0.04,
  ap_coupling  = -0.6,           
  
  # --- Environmental design ---
  n_env        = 12,
  sigma_E      = 3.0,
  env_range    = "wide",         
  
  # --- Index construction ---
  rho          = 0.95,
  kappa        = 1.0,
  
  # --- Residual and testing ---
  sigma_res    = 1.0,
  alpha_bonf   = 0.05,
  hetero_res   = FALSE,          
  
  # --- GWAS method ---
  gwas_method  = "simple",       # "simple", "pca", or "lmm"
  n_pc         = 3,              
  
  # --- Missing data and replication ---
  pct_missing  = 0.0,           
  n_reps       = 1L,            
  
  # --- LD detection window ---
  ld_window    = 5L,            
  
  # --- Nonlinear reaction norm ---
  nonlinear    = FALSE,
  var_a_quad   = 0.0,
  
  # --- Genomic prediction CV toggle ---
  run_gp       = TRUE,
  
  # --- v10 additions. All are off here, so default_params reproduce v9.0.0. ---
  # No-intercept slope QTL class (Section 4B). n_so_qtl = 0 switches it off.
  n_so_qtl     = 0L,
  var_a_so     = 0.04,
  so_axis      = "index",        # "index": effect along x; "driver": along h
  so_center    = "intercept",    # "intercept": zero effect at x = 0;
                                 # "mean": zero effect on the environment-averaged mean
  # Simplified core simulation (Section 2B).
  core_sim     = FALSE,          # independent loci, fixed allele frequency, no LD
  core_maf     = 0.30,
  fixed_effects = FALSE,         # effect magnitude fixed at sqrt(var / 2f(1-f)), random sign
  # Optional per-replicate outputs.
  predict_power = FALSE,         # analytic power prediction (Section 7C)
  pve_cols     = FALSE,          # per-locus variance explained by class (Section 7C)
  audit_cols   = FALSE           # columns used by reaction_norm_audit.R
)

# Standard and Micro-Resolution Rho Grids
rho_grid_fine   <- unique(sort(c(seq(0.00, 0.99, by = 0.01), 
                                 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 
                                 0.993, 0.996, 0.999)))
rho_grid_coarse <- c(0.00, 0.50, 0.70, 0.85, 0.95, 0.99, 0.996)
rho_grid_sensitivity <- c(0.00, 0.50, 0.85, 0.95, 0.996)

# ============================================================================
# SECTION 2: GENOTYPE SIMULATION
# ============================================================================

simulate_genotypes <- function(n, m, maf_min, n_mean, n_shared, n_indep,
                               n_cn = 0, n_ap = 0, pop_structure = "diverse") {
  n_qtl <- n_mean + n_shared + n_indep + n_cn + n_ap
  
  if (pop_structure == "diverse") {
    freq <- runif(m, min = maf_min, max = 0.5)
    block_size <- 20L
    ar1_rho    <- 0.7
    
    ar1_cor <- ar1_rho ^ abs(outer(seq_len(block_size), seq_len(block_size), "-"))
    L_full  <- chol(ar1_cor)   
    
    n_blocks   <- ceiling(m / block_size)
    n_full     <- m %/% block_size         
    n_leftover <- m %% block_size          
    
    make_haplotype <- function() {
      E <- matrix(rnorm(n * m), n, m)      
      if (n_full > 0) {
        for (b in seq_len(n_full)) {
          cols <- ((b - 1L) * block_size + 1L):(b * block_size)
          E[, cols] <- E[, cols] %*% L_full   
        }
      }
      if (n_leftover > 0) {
        cols_last <- (n_full * block_size + 1L):m
        L_part <- chol(ar1_cor[seq_len(n_leftover), seq_len(n_leftover)])
        E[, cols_last] <- E[, cols_last] %*% L_part
      }
      U <- pnorm(E)
      (U < matrix(freq, nrow = n, ncol = m, byrow = TRUE)) + 0L
    }
    
    hap1 <- make_haplotype()
    hap2 <- make_haplotype()
    X <- hap1 + hap2  
    
  } else if (pop_structure == "biparental") {
    freq <- rep(0.5, m)
    parent1 <- rbinom(m, 2, 0.95)  
    parent2 <- rbinom(m, 2, 0.05)  
    X <- matrix(NA_integer_, n, m)
    for (i in seq_len(n)) {
      rec_points1 <- sort(sample(m, rpois(1, m * 0.02)))
      rec_points2 <- sort(sample(m, rpois(1, m * 0.02)))
      gam1 <- parent1; gam2 <- parent2
      
      if (length(rec_points1) > 0) {
        state <- TRUE
        for (rp in rec_points1) {
          if (state) gam1[rp:m] <- parent2[rp:m] else gam1[rp:m] <- parent1[rp:m]
          state <- !state
        }
      }
      if (length(rec_points2) > 0) {
        state <- TRUE
        for (rp in rec_points2) {
          if (state) gam2[rp:m] <- parent1[rp:m] else gam2[rp:m] <- parent2[rp:m]
          state <- !state
        }
      }
      X[i, ] <- pmin(gam1, 1L) + pmin(gam2, 1L)
    }
  } else if (pop_structure == "nam") {
    n_fam  <- 10
    n_per  <- ceiling(n / n_fam)
    freq   <- runif(m, min = maf_min, max = 0.5)
    common <- rbinom(m, 2, 0.5)  
    Xlist  <- vector("list", n_fam)
    for (f in seq_len(n_fam)) {
      diverse_parent <- rbinom(m, 2, freq)
      Xf <- matrix(NA_integer_, n_per, m)
      for (i in seq_len(n_per)) {
        rec <- sort(sample(m, rpois(1, m * 0.02)))
        gam1 <- common; gam2 <- diverse_parent
        if (length(rec) > 0) {
          state <- TRUE
          for (rp in rec) {
            if (state) gam1[rp:m] <- diverse_parent[rp:m] else gam1[rp:m] <- common[rp:m]
            state <- !state
          }
        }
        Xf[i, ] <- pmin(gam1, 1L) + pmin(gam2, 1L)
      }
      Xlist[[f]] <- Xf
    }
    X <- do.call(rbind, Xlist)[seq_len(n), , drop = FALSE]
  }
  
  X_c  <- scale(X, center = TRUE, scale = FALSE)
  freq <- colMeans(X) / 2
  
  pool <- sample(m, n_qtl, replace = FALSE)
  idx  <- list(
    mean   = pool[seq_len(n_mean)],
    shared = pool[n_mean + seq_len(n_shared)],
    indep  = if (n_indep > 0) pool[n_mean + n_shared + seq_len(n_indep)] else integer(0),
    cn     = if (n_cn > 0)    pool[n_mean + n_shared + n_indep + seq_len(n_cn)] else integer(0),
    ap     = if (n_ap > 0)    pool[n_mean + n_shared + n_indep + n_cn + seq_len(n_ap)] else integer(0)
  )
  
  list(X = X_c, freq = freq, qtl_idx = idx, pop_structure = pop_structure)
}

# ============================================================================
# SECTION 2B: SIMPLIFIED CORE GENOTYPES (v10; used only when core_sim = TRUE)
# ============================================================================
# Independent loci (no LD), one fixed allele frequency for every marker. The QTL
# pool is drawn exactly as in simulate_genotypes(), so the class layout is the same.

simulate_genotypes_core <- function(n, m, maf, n_mean, n_shared, n_indep,
                                    n_cn = 0, n_ap = 0) {
  n_qtl <- n_mean + n_shared + n_indep + n_cn + n_ap
  X <- matrix(rbinom(n * m, 2L, maf), n, m)
  X_c  <- scale(X, center = TRUE, scale = FALSE)
  freq <- colMeans(X) / 2
  pool <- sample(m, n_qtl, replace = FALSE)
  idx  <- list(
    mean   = pool[seq_len(n_mean)],
    shared = pool[n_mean + seq_len(n_shared)],
    indep  = if (n_indep > 0) pool[n_mean + n_shared + seq_len(n_indep)] else integer(0),
    cn     = if (n_cn > 0)    pool[n_mean + n_shared + n_indep + seq_len(n_cn)] else integer(0),
    ap     = if (n_ap > 0)    pool[n_mean + n_shared + n_indep + n_cn + seq_len(n_ap)] else integer(0)
  )
  list(X = X_c, freq = freq, qtl_idx = idx, pop_structure = "core")
}

# ============================================================================
# SECTION 3: PHENOTYPE GENERATION
# ============================================================================

generate_phenotypes <- function(geno, params) {
  n   <- params$n_geno
  m_e <- params$n_env
  X   <- geno$X
  idx <- geno$qtl_idx
  
  sfx <- function(nq, vq, fq) {
    if (nq == 0) return(numeric(0))
    sdv <- sqrt(vq / pmax(2 * fq * (1 - fq), 0.01))
    # v10: fixed effect magnitudes with random sign (simplified core simulation).
    if (isTRUE(params$fixed_effects)) return(sample(c(-1, 1), nq, replace = TRUE) * sdv)
    rnorm(nq, 0, sdv)
  }
  
  a_mean   <- sfx(params$n_mean_qtl, params$var_a_mean, geno$freq[idx$mean])
  a_sh_int <- sfx(params$n_shared_qtl, params$var_a_shared, geno$freq[idx$shared])
  a_sh_slp <- params$ds_coupling * a_sh_int
  a_indep  <- sfx(params$n_indep_qtl, params$var_a_indep, if (length(idx$indep) > 0) geno$freq[idx$indep] else numeric(0))
  a_cn     <- sfx(params$n_cn_qtl, params$var_a_cn, if (length(idx$cn) > 0) geno$freq[idx$cn] else numeric(0))
  a_ap_int <- sfx(params$n_ap_qtl, params$var_a_ap, if (length(idx$ap) > 0) geno$freq[idx$ap] else numeric(0))
  a_ap_slp <- params$ap_coupling * a_ap_int 
  
  a_quad <- if (params$nonlinear && params$var_a_quad > 0) {
    sfx(params$n_mean_qtl, params$var_a_quad, geno$freq[idx$mean])
  } else rep(0, params$n_mean_qtl)
  
  g_mean   <- as.vector(X[, idx$mean, drop = FALSE] %*% a_mean)
  g_sh_int <- as.vector(X[, idx$shared, drop = FALSE] %*% a_sh_int)
  g_sh_slp <- as.vector(X[, idx$shared, drop = FALSE] %*% a_sh_slp)
  g_indep  <- if (params$n_indep_qtl > 0) as.vector(X[, idx$indep, drop = FALSE] %*% a_indep) else rep(0, n)
  g_cn     <- if (params$n_cn_qtl > 0) as.vector(X[, idx$cn, drop = FALSE] %*% a_cn) else rep(0, n)
  g_ap_int <- if (params$n_ap_qtl > 0) as.vector(X[, idx$ap, drop = FALSE] %*% a_ap_int) else rep(0, n)
  g_ap_slp <- if (params$n_ap_qtl > 0) as.vector(X[, idx$ap, drop = FALSE] %*% a_ap_slp) else rep(0, n)
  g_quad   <- if (params$nonlinear && params$var_a_quad > 0) as.vector(X[, idx$mean, drop = FALSE] %*% a_quad) else rep(0, n)
  
  true_int <- g_mean + g_sh_int + g_ap_int
  
  h_sd <- if (params$env_range == "narrow") 0.3 else 1.0
  h    <- rnorm(m_e, 0, h_sd)
  
  z_raw <- rnorm(m_e, 0, 1.0)
  z     <- residuals(lm(z_raw ~ h))
  z     <- z / sd(z)
  
  e_j <- params$sigma_E * h
  cn_indicator <- as.numeric(h > params$cn_threshold)
  
  mu <- 50
  Y  <- matrix(NA_real_, n, m_e)
  
  if (isTRUE(params$hetero_res)) {
    shape <- 4; rate <- shape / params$sigma_res
    sig_j <- sqrt(rgamma(m_e, shape = shape, rate = rate))
  } else {
    sig_j <- rep(params$sigma_res, m_e)
  }
  
  for (j in seq_len(m_e)) {
    Y[, j] <- mu + true_int + e_j[j] + g_sh_slp * h[j] + g_indep * z[j] +
      g_cn * cn_indicator[j] + g_ap_slp * h[j] + g_quad * h[j]^2 + rnorm(n, 0, sig_j[j])
  }
  
  n_reps_param <- if (!is.null(params$n_reps)) params$n_reps else 1L
  if (n_reps_param > 1L) {
    for (rep_k in 2:n_reps_param) {
      Y_rep <- matrix(NA_real_, n, m_e)
      for (j in seq_len(m_e)) {
        Y_rep[, j] <- mu + true_int + e_j[j] + g_sh_slp * h[j] + g_indep * z[j] +
          g_cn * cn_indicator[j] + g_ap_slp * h[j] + g_quad * h[j]^2 + rnorm(n, 0, sig_j[j])
      }
      Y <- Y + Y_rep
    }
    Y <- Y / n_reps_param 
  }
  
  pct_miss <- if (!is.null(params$pct_missing)) params$pct_missing else 0
  if (pct_miss > 0) {
    n_miss <- round(n * m_e * pct_miss)
    if (n_miss > 0) Y[sample(n * m_e, n_miss, replace = FALSE)] <- NA_real_
    
    for (i in seq_len(n)) {
      n_obs <- sum(!is.na(Y[i, ]))
      if (n_obs < 3) {
        na_cols <- which(is.na(Y[i, ]))
        restore <- na_cols[seq_len(min(3 - n_obs, length(na_cols)))]
        for (rc in restore) {
          Y[i, rc] <- mu + true_int[i] + e_j[rc] + g_sh_slp[i] * h[rc] + g_indep[i] * z[rc] +
            g_cn[i] * cn_indicator[rc] + g_ap_slp[i] * h[rc] + g_quad[i] * h[rc]^2 + rnorm(1, 0, sig_j[rc])
        }
      }
    }
  }
  
  list(Y = Y, env_mean = colMeans(Y, na.rm = TRUE), h = h, z = z, e_j = e_j,
       cn_indicator = cn_indicator, true_intercept = true_int, sig_j = sig_j,
       true_slope_h = g_sh_slp + g_ap_slp, true_slope_z = g_indep,
       true_cn = g_cn, true_quad = g_quad,
       effects = list(a_mean = a_mean, a_sh_int = a_sh_int, a_sh_slp = a_sh_slp,
                      a_indep = a_indep, a_cn = a_cn, a_ap_int = a_ap_int, 
                      a_ap_slp = a_ap_slp, a_quad = a_quad))
}

# ============================================================================
# SECTION 4 & 5: INDEX & FINLAY-WILKINSON
# ============================================================================

construct_index <- function(pheno, rho, kappa = 1.0) {
  em_std <- as.vector(scale(pheno$env_mean))
  z_raw  <- as.vector(scale(pheno$z))
  m_e    <- length(em_std)
  
  z_orth <- z_raw - as.numeric(crossprod(z_raw, em_std) / crossprod(em_std)) * em_std
  zn <- sqrt(sum(z_orth^2) / (m_e - 1))
  if (zn > 1e-10) z_orth <- z_orth / zn else {
    z_orth <- scale(rnorm(m_e))[,1]
  }
  
  w <- rnorm(m_e)
  w <- w - as.numeric(crossprod(w, em_std) / crossprod(em_std)) * em_std
  w <- w - as.numeric(crossprod(w, z_orth) / crossprod(z_orth)) * z_orth
  wn <- sqrt(sum(w^2) / (m_e - 1))
  if (wn > 1e-10) w <- w / wn
  
  tau      <- kappa * sqrt(max(1 - rho^2, 0))
  noise_sd <- sqrt(max(1 - rho^2 - tau^2, 0))
  
  x <- rho * em_std + tau * z_orth + noise_sd * w
  list(x = x, target_rho = rho, actual_rho = cor(x, pheno$env_mean), tau = tau)
}

# ============================================================================
# SECTION 4B: NO-INTERCEPT SLOPE QTL CLASS (v10; active only when n_so_qtl > 0)
# ============================================================================
# A slope-only locus k adds g_k,i * v_j to genotype i in environment j, where v is
# an environmental axis and g_k,i = X_ik * a_k. Two axes and two centerings:
#   so_axis = "index"  : v is built from the index x itself (the CERIS regressor).
#   so_axis = "driver" : v is built from the environmental driver h.
#   so_center = "intercept" : v is shifted so that its least-squares regression
#       on x passes through the origin. The locus then changes the slope on x and
#       has exactly zero effect on the intercept at x = 0. For so_axis = "index"
#       this is v = x.
#   so_center = "mean" : v is centered at its mean over the environments, so the
#       locus has exactly zero effect on the environment-averaged genotype mean.
#   The two centerings coincide when mean(x) = 0. In this simulation mean(x) = 0
#   when kappa = 1 (the default), because x = rho*y + tau*z_perp + nu*w, y and
#   z_perp are centered, and nu = 0 (up to floating-point rounding, which leaves
#   |mean(x)| below about 1e-8). For kappa < 1, mean(x) = nu*mean(w).
# Because genotypes are column-centered, sum_i g_k,i = 0, so the class leaves the
# environment means, and therefore the index x, unchanged. The term is added after
# the index is built, which is why so_axis = "index" is well defined.
# Positions and effects are drawn from an independent RNG substream
# (seed + .SO_SEED_OFFSET) and the main stream is restored afterwards, so a
# replicate with the class switched on is identical to the same replicate with it
# off in every other locus, effect, environment, residual and index draw.

.SO_SEED_OFFSET <- 700000000L

with_substream <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}

so_axis_vector <- function(axis, center, pheno, idx) {
  x  <- idx$x
  v  <- if (identical(as.character(axis), "driver")) pheno$h else x
  xc <- x - mean(x)
  b_vx <- sum(xc * v) / sum(xc^2)
  if (identical(as.character(center), "mean")) v - mean(v) else v - mean(v) + b_vx * mean(x)
}

add_slope_only_qtl <- function(geno, pheno, idx, params, seed) {
  n_so <- if (is.null(params$n_so_qtl)) 0L else as.integer(params$n_so_qtl)
  if (is.na(n_so) || n_so <= 0L) {
    geno$qtl_idx$so <- integer(0)
    return(list(geno = geno, pheno = pheno))
  }
  avail <- setdiff(seq_len(ncol(geno$X)), unlist(geno$qtl_idx))
  drawn <- with_substream(seed + .SO_SEED_OFFSET, {
    pos <- avail[sample.int(length(avail), n_so)]
    f   <- geno$freq[pos]
    sdv <- sqrt(params$var_a_so / pmax(2 * f * (1 - f), 0.01))
    a   <- if (isTRUE(params$fixed_effects)) sample(c(-1, 1), n_so, replace = TRUE) * sdv else rnorm(n_so, 0, sdv)
    list(pos = pos, a = a)
  })
  v    <- so_axis_vector(params$so_axis, params$so_center, pheno, idx)
  g_so <- as.vector(geno$X[, drawn$pos, drop = FALSE] %*% drawn$a)
  pheno$Y <- pheno$Y + outer(g_so, v)
  pheno$g_so <- g_so; pheno$so_vec <- v; pheno$effects$a_so <- drawn$a
  geno$qtl_idx$so <- drawn$pos
  list(geno = geno, pheno = pheno)
}

# --- Section 5 (continued): reaction norm estimation on the index x ---
# Slope b_i = sum_j (x_j - mean(x)) Y_ij / sum_j (x_j - mean(x))^2 and intercept
# a_i = mean_j(Y_ij) - b_i * mean(x), so the intercept is the fitted value at x = 0.
# It equals the environment-averaged genotype mean only when mean(x) = 0.

estimate_reaction_norms <- function(Y, x) {
  n <- nrow(Y); m_e <- ncol(Y)
  xc <- x - mean(x)
  slopes <- numeric(n); intercepts <- numeric(n); mean_perfs <- numeric(n)
  
  if (!anyNA(Y)) {
    vx <- sum(xc^2)
    ym <- rowMeans(Y)
    slopes <- as.vector((Y - ym) %*% xc) / vx
    intercepts <- ym - slopes * mean(x)
    mean_perfs <- ym
  } else {
    mx <- mean(x)
    for (i in seq_len(n)) {
      obs <- !is.na(Y[i, ])
      if (sum(obs) < 2) {
        slopes[i] <- NA_real_; intercepts[i] <- NA_real_
        mean_perfs[i] <- if (any(obs)) mean(Y[i, obs]) else NA_real_
        next
      }
      yi <- Y[i, obs]; xi <- x[obs]; xci <- xi - mean(xi)
      slopes[i] <- sum(yi * xci) / sum(xci^2)
      intercepts[i] <- mean(yi) - slopes[i] * mean(xi)
      mean_perfs[i] <- mean(yi)
    }
  }
  data.frame(geno_id = seq_len(n), intercept = intercepts, slope = slopes, mean_perf = mean_perfs)
}

# ============================================================================
# SECTION 6: FAST MARGINAL GWAS (Simple, PCA, and LMM)
# ============================================================================

fast_gwas <- function(X, y) {
  n <- length(y)
  yc <- y - mean(y)
  ssx <- colSums(X^2)
  ssy <- sum(yc^2)
  xy <- as.vector(crossprod(X, yc))
  
  mono <- ssx < .Machine$double.eps
  ssx[mono] <- 1  
  
  beta <- xy / ssx
  ssr  <- pmax(ssy - beta * xy, .Machine$double.eps)
  se   <- sqrt((ssr / (n - 2)) / ssx)
  tt   <- beta / se
  pv   <- 2 * pt(abs(tt), df = n - 2, lower.tail = FALSE)
  
  beta[mono] <- 0; se[mono] <- Inf; tt[mono] <- 0; pv[mono] <- 1
  fdr <- p.adjust(pv, method = "BH")
  
  data.frame(marker = seq_along(beta), beta = beta, se = se, t_stat = tt, p_value = pv, fdr = fdr)
}

fast_gwas_pca <- function(X, y, n_pc = 3, PCs = NULL) {
  n <- length(y)
  if (n_pc > 0 && n > n_pc + 2 && !is.null(PCs)) {
    Q <- qr(cbind(1, PCs))
    y_r <- qr.resid(Q, y)
    X_r <- qr.resid(Q, X)
    df_res <- n - n_pc - 2
  } else {
    y_r <- y - mean(y)
    X_r <- X
    df_res <- n - 2
  }
  
  ssx <- colSums(X_r^2)
  ssy <- sum(y_r^2)
  xy  <- as.vector(crossprod(X_r, y_r))
  
  mono <- ssx < .Machine$double.eps
  ssx[mono] <- 1
  
  beta <- xy / ssx
  ssr  <- pmax(ssy - beta * xy, .Machine$double.eps)
  se   <- sqrt((ssr / max(df_res, 1)) / ssx)
  tt   <- beta / se
  pv   <- 2 * pt(abs(tt), df = max(df_res, 1), lower.tail = FALSE)
  
  beta[mono] <- 0; se[mono] <- Inf; tt[mono] <- 0; pv[mono] <- 1
  fdr  <- p.adjust(pv, method = "BH")
  
  data.frame(marker = seq_along(beta), beta = beta, se = se, t_stat = tt, p_value = pv, fdr = fdr)
}

# LMM GWAS: rotation by the inverse Cholesky factor L^{-1}
fast_gwas_lmm <- function(X, y, K = NULL, lambda = 1.0) {
  n <- length(y)
  if (is.null(K)) K <- tcrossprod(X) / ncol(X)
  
  V <- K
  diag(V) <- diag(V) + lambda
  
  U <- tryCatch(chol(V), error = function(e) NULL)
  
  if (is.null(U)) {
    return(fast_gwas(X, y))
  }
  
  # Correct rotation: L^{-1}y where L = t(U) (lower triangular)
  # forwardsolve(L, y) computes L^{-1}y
  y_c <- y - mean(y)
  L <- t(U)
  y_star <- forwardsolve(L, y_c)
  X_star <- forwardsolve(L, X)
  
  ssx <- colSums(X_star^2)
  ssy <- sum(y_star^2)
  xy  <- as.vector(crossprod(X_star, y_star))
  
  mono <- ssx < .Machine$double.eps
  ssx[mono] <- 1
  
  beta <- xy / ssx
  ssr  <- pmax(ssy - beta * xy, .Machine$double.eps)
  se   <- sqrt((ssr / (n - 2)) / ssx)
  tt   <- beta / se
  pv   <- 2 * pt(abs(tt), df = n - 2, lower.tail = FALSE)
  
  beta[mono] <- 0; se[mono] <- Inf; tt[mono] <- 0; pv[mono] <- 1
  fdr  <- p.adjust(pv, method = "BH")
  
  data.frame(marker = seq_along(beta), beta = beta, se = se, t_stat = tt, p_value = pv, fdr = fdr)
}

# ============================================================================
# SECTION 6B: GENOMIC PREDICTION (Optimized Ridge Regression)
# ============================================================================

genomic_prediction_cv <- function(X, y_int, y_slope, n_folds = 5) {
  n <- nrow(X); m <- ncol(X)
  if (n > 2500) return(list(gp_acc_int = NA, gp_acc_slope = NA, gp_acc_int2 = NA, gp_acc_slope2 = NA))
  
  folds <- sample(rep(seq_len(n_folds), length.out = n))
  pred_int <- numeric(n); pred_slope <- numeric(n)
  
  K <- tcrossprod(X) / m
  lambda <- n / m   
  
  for (f in seq_len(n_folds)) {
    test  <- which(folds == f); train <- which(folds != f)
    n_tr  <- length(train)
    
    K_tr <- K[train, train]
    K_tr[matrix(c(1:n_tr, 1:n_tr), ncol = 2)] <- K_tr[matrix(c(1:n_tr, 1:n_tr), ncol = 2)] + lambda
    
    U <- tryCatch(chol(K_tr), error = function(e) NULL)
    if (is.null(U)) next 
    
    y_c_int <- y_int[train] - mean(y_int[train])
    alpha_int <- backsolve(U, forwardsolve(t(U), y_c_int))
    pred_int[test] <- as.vector(K[test, train] %*% alpha_int) + mean(y_int[train])
    
    y_c_slp <- y_slope[train] - mean(y_slope[train])
    alpha_slp <- backsolve(U, forwardsolve(t(U), y_c_slp))
    pred_slope[test] <- as.vector(K[test, train] %*% alpha_slp) + mean(y_slope[train])
  }
  
  list(gp_acc_int   = cor(pred_int, y_int, use = "complete.obs"),
       gp_acc_slope = cor(pred_slope, y_slope, use = "complete.obs"),
       gp_acc_int2  = cor(pred_int, y_int, use = "complete.obs")^2,
       gp_acc_slope2 = cor(pred_slope, y_slope, use = "complete.obs")^2)
}

# ============================================================================
# SECTION 6C: AMMI-1 DECOMPOSITION
# ============================================================================

estimate_ammi <- function(Y) {
  n <- nrow(Y); m_e <- ncol(Y)
  grand_mu <- mean(Y, na.rm = TRUE)
  g_i <- rowMeans(Y, na.rm = TRUE) - grand_mu
  e_j <- colMeans(Y, na.rm = TRUE) - grand_mu
  
  GxE <- Y - outer(g_i, rep(1, m_e)) - rep(e_j, each = n) - grand_mu
  if (anyNA(GxE)) GxE[is.na(GxE)] <- 0
  
  sv <- svd(GxE, nu = 1, nv = 1)
  gamma <- sv$u[, 1] * sqrt(sv$d[1])
  ss_gxe <- sum(GxE^2, na.rm = TRUE)
  pct_ipca1 <- if (ss_gxe > 0) sv$d[1]^2 / ss_gxe * 100 else NA_real_
  
  data.frame(geno_id = seq_len(n), intercept = g_i + grand_mu, ipca1 = gamma, 
             mean_perf = rowMeans(Y, na.rm = TRUE), pct_ipca1 = pct_ipca1)
}

# ============================================================================
# SECTION 7: METRICS COMPILER
# ============================================================================

compute_metrics <- function(gi, gs, rn, qi, nm, pheno_Y, ab = 0.05, ld_window = 5L) {
  th <- ab / nm; fdr_th <- 0.05
  
  si <- na.omit(gi$marker[gi$p_value < th])
  ss <- na.omit(gs$marker[gs$p_value < th])
  si_fdr <- na.omit(if ("fdr" %in% names(gi)) gi$marker[gi$fdr < fdr_th] else si)
  ss_fdr <- na.omit(if ("fdr" %in% names(gs)) gs$marker[gs$fdr < fdr_th] else ss)
  
  pw_window <- function(qtl_pos, sig_set) {
    if (length(qtl_pos) == 0 || length(sig_set) == 0) return(NA_real_)
    mean(vapply(qtl_pos, function(q) any(abs(sig_set - q) <= ld_window), logical(1)))
  }
  pw <- function(q, s) if (length(q) > 0) mean(q %in% s) else NA_real_
  
  nss <- length(ss); nsi <- length(si)
  ov  <- length(intersect(ss, si)); us  <- length(setdiff(ss, si))
  
  all_qtl <- sort(unique(c(qi$shared, qi$indep, qi$cn, qi$ap, qi$so)))
  window_det_int <- if (length(all_qtl) > 0) vapply(all_qtl, function(q) any(abs(si - q) <= ld_window), logical(1)) else logical(0)
  window_det_slp <- if (length(all_qtl) > 0) vapply(all_qtl, function(q) any(abs(ss - q) <= ld_window), logical(1)) else logical(0)
  n_win_both <- sum(window_det_int & window_det_slp); n_win_slp <- sum(window_det_slp)
  
  null_markers <- setdiff(seq_len(nm), unlist(qi))
  chi2_med <- qchisq(0.5, df = 1)
  lambda_int <- if (length(null_markers) > 10) median(gi$t_stat[null_markers]^2, na.rm = TRUE) / chi2_med else NA_real_
  lambda_slope <- if (length(null_markers) > 10) median(gs$t_stat[null_markers]^2, na.rm = TRUE) / chi2_med else NA_real_
  
  csi <- cor(rn$slope, rn$intercept, use = "complete.obs")
  csm <- cor(rn$slope, rn$mean_perf, use = "complete.obs")
  
  rank_cors <- if (ncol(pheno_Y) >= 2) apply(combn(ncol(pheno_Y), 2), 2, function(ij) cor(pheno_Y[, ij[1]], pheno_Y[, ij[2]], method = "spearman")) else numeric(0)
  
  # FDR-based power metrics
  fdr_pw_shared <- pw(qi$shared, ss_fdr)
  fdr_pw_indep  <- pw(qi$indep, ss_fdr)
  
  list(
    n_sig_slope = nss, n_sig_int = nsi, n_overlap = ov, n_unique_slope = us,
    overlap_prop = if (nss > 0) ov / nss else NA_real_,
    unique_prop = if (nss > 0) us / nss else NA_real_,
    cor_slope_int = csi, cor_slope_mean = csm,
    cor_logp = cor(-log10(gi$p_value), -log10(gs$p_value), method = "spearman"),
    power_mean_in_int = pw(qi$mean, si), power_shared_in_int = pw(qi$shared, si),
    power_shared_in_slope = pw(qi$shared, ss), power_indep_in_slope = pw(qi$indep, ss),
    power_cn_in_slope = pw(qi$cn, ss), power_ap_in_slope = pw(qi$ap, ss),
    fpr_int = pw(null_markers, si), fpr_slope = pw(null_markers, ss),
    lambda_int = lambda_int, lambda_slope = lambda_slope,
    fdr_overlap_prop = if (length(ss_fdr) > 0) length(intersect(ss_fdr, si_fdr)) / length(ss_fdr) else NA_real_,
    fdr_unique_prop = if (length(ss_fdr) > 0) length(setdiff(ss_fdr, si_fdr)) / length(ss_fdr) else NA_real_,
    fdr_power_shared = fdr_pw_shared,
    fdr_power_indep = fdr_pw_indep,
    median_rank_cor = if (length(rank_cors) > 0) median(rank_cors) else NA_real_,
    prop_low_rank_cor = if (length(rank_cors) > 0) mean(rank_cors < 0.8) else NA_real_,
    win_power_shared = pw_window(qi$shared, ss), win_power_indep = pw_window(qi$indep, ss),
    win_overlap_prop = if (n_win_slp > 0) n_win_both / n_win_slp else NA_real_,
    win_unique_prop = if (n_win_slp > 0) 1 - (n_win_both / n_win_slp) else NA_real_
  )
}

# ============================================================================
# SECTION 7B: QUANTITATIVE GENETICS THEORY METRICS
# ============================================================================

compute_qg_metrics <- function(pheno, idx, rn, geno, params) {
  valid <- !is.na(rn$slope) & !is.na(rn$intercept)
  if (sum(valid) < 10) return(list(r_predicted = NA, r_true_genetic = NA, V_b_cond_a = NA, pct_V_b_unique = NA, mean_type_b = NA, min_type_b = NA, r_g_via_lande = NA, n_eff_grn = NA, h2_b_narrow = NA, R_b_std = NA, R_b_independent = NA, R_b_ind_std = NA, n_eff_env = NA, pct_pc1 = NA, noise_inflation = NA))
  
  true_a <- pheno$true_intercept; true_s <- pheno$true_slope_h
  true_t <- pheno$true_slope_z; true_cn <- pheno$true_cn
  V_a <- var(true_a); V_s <- var(true_s); V_t <- var(true_t); C_as <- cov(true_a, true_s)
  rho <- idx$actual_rho; tau <- idx$tau
  
  V_bx_true <- rho^2 * V_s + tau^2 * V_t; C_abx <- rho * C_as
  r_predicted <- if (V_a > 1e-12 && V_bx_true > 1e-12) C_abx / sqrt(V_a * V_bx_true) else NA_real_
  
  true_bx <- rho * true_s + tau * true_t
  r_true_genetic <- if (var(true_bx) > 1e-12) cor(true_a, true_bx) else NA_real_
  
  V_b_est <- var(rn$slope, na.rm = TRUE); r_ab_est <- cor(rn$slope, rn$intercept, use = "complete.obs")
  V_b_cond_a <- V_b_est * (1 - r_ab_est^2)
  pct_V_b_unique <- if (V_b_est > 1e-12) V_b_cond_a / V_b_est * 100 else NA_real_
  
  Y_adj <- pheno$Y - rep(colMeans(pheno$Y, na.rm = TRUE), each = nrow(pheno$Y))
  type_b_cors <- apply(combn(ncol(pheno$Y), 2), 2, function(ij) cor(Y_adj[, ij[1]], Y_adj[, ij[2]], use = "complete.obs"))
  mean_type_b <- mean(type_b_cors)
  min_type_b  <- min(type_b_cors)
  
  G_env <- cov(Y_adj, use = "pairwise.complete.obs")
  x_c <- idx$x - mean(idx$x); W <- cbind(rep(1, length(idx$x)), x_c)
  WtW_inv <- tryCatch(solve(crossprod(W)), error = function(e) NULL)
  G_rn_vl <- if (!is.null(WtW_inv)) WtW_inv %*% t(W) %*% G_env %*% W %*% WtW_inv else NULL
  r_g_vl <- if (!is.null(G_rn_vl) && G_rn_vl[1, 1] > 1e-12 && G_rn_vl[2, 2] > 1e-12) G_rn_vl[1, 2] / sqrt(G_rn_vl[1, 1] * G_rn_vl[2, 2]) else NA_real_
  eig_grn <- if (!is.null(G_rn_vl)) pmax(eigen(G_rn_vl, symmetric = TRUE, only.values = TRUE)$values, 0) else c(NA, NA)
  n_eff_grn <- if (!any(is.na(eig_grn)) && sum(eig_grn) > 1e-12) (sum(eig_grn))^2 / sum(eig_grn^2) else NA_real_
  
  sigma_b <- sd(rn$slope, na.rm = TRUE)
  h2_b_narrow <- if (var(true_bx) > 1e-12 && V_b_est > 1e-12) min(var(true_bx) / V_b_est, 1.0) else NA_real_
  R_b <- if (!is.na(h2_b_narrow)) h2_b_narrow * 1.755 * sigma_b else NA_real_
  R_b_std <- if (!is.na(R_b) && sigma_b > 1e-12) R_b / sigma_b else NA_real_
  
  R_b_independent <- if (!is.na(R_b) && !is.na(r_ab_est)) R_b * sqrt(max(1 - r_ab_est^2, 0)) else NA_real_
  R_b_ind_std <- if (!is.na(R_b_independent) && sigma_b > 1e-12) R_b_independent / sigma_b else NA_real_
  
  eig_genv <- pmax(eigen(G_env, symmetric = TRUE, only.values = TRUE)$values, 0)
  n_eff_env <- if (sum(eig_genv) > 1e-12) (sum(eig_genv))^2 / sum(eig_genv^2) else NA_real_
  pct_pc1 <- if (sum(eig_genv) > 1e-12) eig_genv[1] / sum(eig_genv) * 100 else NA_real_
  
  list(r_predicted = r_predicted, r_true_genetic = r_true_genetic, V_b_cond_a = V_b_cond_a, pct_V_b_unique = pct_V_b_unique, mean_type_b = mean_type_b, min_type_b = min_type_b, r_g_via_lande = r_g_vl, n_eff_grn = n_eff_grn, h2_b_narrow = h2_b_narrow, R_b_std = R_b_std, R_b_independent = R_b_independent, R_b_ind_std = R_b_ind_std, n_eff_env = n_eff_env, pct_pc1 = pct_pc1, noise_inflation = if (!is.na(r_predicted) && abs(r_predicted) > 1e-6) r_ab_est / r_predicted else NA_real_)
}

# ============================================================================
# SECTION 7C: ANALYTIC POWER, VARIANCE EXPLAINED AND v10 OPTIONAL COLUMNS
# ============================================================================
# See Appendix A of the paper. The slope and intercept are linear contrasts of each
# genotype's environment-specific values: b_i = sum_j c_j Ybar_ij and
# a_i = sum_j d_j Ybar_ij, with c_j = (x_j - mean(x)) / S_xx and
# d_j = 1/m - mean(x) * c_j. A locus k with per-allele effect A_kj in environment j
# has per-allele slope effect e_k = sum_j c_j A_kj and intercept effect
# f_k = sum_j d_j A_kj. With genotype variance v_k = Var(X_k), the marginal
# t-test of slope on X_k has squared non-centrality
#   ncp2_k = n * v_k e_k^2 / (D_b - v_k e_k^2),
#   D_b = Var_i(sum_j c_j G_ij) + sum_j c_j^2 sigma_j^2 / n_rep,
# where G_ij is the true genotypic value. Power is P(|T_{n-2}(sqrt(ncp2))| > t*)
# with t* the two-sided Bonferroni critical value at alpha_bonf / n_markers.
# Nothing here is fitted to the GWAS outcome.

locus_effect_matrix <- function(geno, pheno, params) {
  e <- pheno$effects; qi <- geno$qtl_idx; m_e <- length(pheno$h)
  h <- pheno$h; z <- pheno$z; cn <- pheno$cn_indicator
  one <- rep(1, m_e)
  blocks <- list(
    mean   = list(pos = qi$mean,   A = outer(e$a_mean, one) + outer(e$a_quad, h^2)),
    shared = list(pos = qi$shared, A = outer(e$a_sh_int, one) + outer(e$a_sh_slp, h)),
    indep  = list(pos = qi$indep,  A = outer(e$a_indep, z)),
    cn     = list(pos = qi$cn,     A = outer(e$a_cn, cn)),
    ap     = list(pos = qi$ap,     A = outer(e$a_ap_int, one) + outer(e$a_ap_slp, h))
  )
  if (length(qi$so) > 0) blocks$so <- list(pos = qi$so, A = outer(e$a_so, pheno$so_vec))
  blocks <- blocks[vapply(blocks, function(b) length(b$pos) > 0, logical(1))]
  list(pos = unlist(lapply(blocks, `[[`, "pos"), use.names = FALSE),
       cls = rep(names(blocks), vapply(blocks, function(b) length(b$pos), integer(1))),
       A   = do.call(rbind, lapply(blocks, function(b) matrix(b$A, nrow = length(b$pos)))))
}

analytic_power_conditional <- function(geno, pheno, idx, params, eff = NULL) {
  if (is.null(eff)) eff <- locus_effect_matrix(geno, pheno, params)
  x <- idx$x; m_e <- length(x); n <- nrow(geno$X)
  xc <- x - mean(x); cvec <- xc / sum(xc^2); dvec <- rep(1 / m_e, m_e) - mean(x) * cvec
  Xq <- geno$X[, eff$pos, drop = FALSE]
  vk <- colMeans(Xq^2)
  es <- as.vector(eff$A %*% cvec); ei <- as.vector(eff$A %*% dvec)
  G  <- Xq %*% eff$A
  vpop <- function(u) mean((u - mean(u))^2)
  nr <- if (is.null(params$n_reps)) 1L else params$n_reps
  sig2 <- (if (!is.null(pheno$sig_j)) pheno$sig_j^2 else rep(params$sigma_res^2, m_e)) / nr
  D_b <- vpop(as.vector(G %*% cvec)) + sum(cvec^2 * sig2)
  D_a <- vpop(as.vector(G %*% dvec)) + sum(dvec^2 * sig2)
  ncp2_b <- n * vk * es^2 / pmax(D_b - vk * es^2, 1e-12)
  ncp2_a <- n * vk * ei^2 / pmax(D_a - vk * ei^2, 1e-12)
  df <- n - 2; tcrit <- qt(1 - params$alpha_bonf / (2 * params$n_markers), df)
  pw <- function(l) pt(-tcrit, df, ncp = sqrt(l)) + pt(tcrit, df, ncp = sqrt(l), lower.tail = FALSE)
  list(cls = eff$cls, pos = eff$pos, vk = vk, e_slope = es, e_int = ei,
       ncp2_slope = ncp2_b, ncp2_int = ncp2_a,
       power_slope = pw(ncp2_b), power_int = pw(ncp2_a),
       D_slope = D_b, D_int = D_a, mean_x = mean(x), sum_c2 = sum(cvec^2),
       G = G, cvec = cvec, dvec = dvec)
}

# Per-locus variance explained, averaged within class: in each environment
# (v_k A_kj^2 / Var_i(Y_ij), averaged over j), on the genotype means over
# environments (v_k mean_j(A_kj)^2 / Var_i(mean_j Y_ij)), on the estimated slope
# (v_k e_k^2 / Var(b)) and on the estimated intercept (v_k f_k^2 / Var(a)).
locus_pve <- function(ap, pheno, rn) {
  vpop <- function(u) { u <- u[!is.na(u)]; mean((u - mean(u))^2) }
  Y <- pheno$Y
  vp_env <- apply(Y, 2, vpop)
  vp_gm  <- vpop(rowMeans(Y, na.rm = TRUE))
  eff_env <- sweep(matrix(ap$vk, nrow = length(ap$vk), ncol = length(vp_env)) *
                     (ap$A_full^2), 2, vp_env, "/")
  data.frame(cls = ap$cls,
             pve_env   = rowMeans(eff_env),
             pve_gmean = ap$vk * rowMeans(ap$A_full)^2 / vp_gm,
             pve_slope = ap$vk * ap$e_slope^2 / vpop(rn$slope),
             pve_int   = ap$vk * ap$e_int^2 / vpop(rn$intercept),
             stringsAsFactors = FALSE)
}

# Closed-form expected power (Appendix A of the paper). Power is
# averaged over the sampling distribution of S_h^2, the variance of the driver h
# across the m environments (average_env = TRUE); D and ncp2 are reported at
# E[S_h^2] = sigma_h^2. Other assumptions: the index
# residual is fully informative up to kappa (tau = kappa*sqrt(1 - rho^2)),
# mean(x) = 0, S_xx = m - 1, and loci are unlinked. Homoscedastic residuals and
# a linear reaction norm are assumed. With random_effects = TRUE each locus
# variance is sigma_a^2 * chi-square(1), as in sfx(); otherwise it is sigma_a^2.
closed_form_power <- function(params, rho = params$rho, kappa = params$kappa,
                              random_effects = !isTRUE(params$fixed_effects),
                              average_env = TRUE, n_env_quad = 41L) {
  n <- params$n_geno; m <- params$n_env; nr <- if (is.null(params$n_reps)) 1L else params$n_reps
  s_h <- if (identical(as.character(params$env_range), "narrow")) 0.3 else 1.0
  tau <- kappa * sqrt(max(1 - rho^2, 0))
  th  <- params$cn_threshold / s_h
  p_cn <- 1 - pnorm(th); phi_cn <- dnorm(th)
  so_axis <- as.character(params$so_axis); n_so <- if (is.null(params$n_so_qtl)) 0L else params$n_so_qtl
  df <- n - 2; tcrit <- qt(1 - params$alpha_bonf / (2 * params$n_markers), df)
  pw <- function(l) pt(-tcrit, df, ncp = sqrt(l)) + pt(tcrit, df, ncp = sqrt(l), lower.tail = FALSE)
  zq <- seq(-8, 8, length.out = 801); wq <- dnorm(zq); wq <- wq / sum(wq)
  epw <- function(v, e2, D) {
    if (v * e2 <= 0) return(pw(0))
    if (!random_effects) return(pw(n * v * e2 / max(D - v * e2, 1e-12)))
    u <- zq^2
    sum(wq * pw(n * v * u * e2 / pmax(D - v * u * e2, 1e-12)))
  }
  # per class: count, sigma_a^2, slope effect^2 and intercept effect^2 per unit allele
  # effect, given the sample variance S2 of the driver h over the m environments
  classes <- function(S2) data.frame(
    cls  = c("mean", "shared", "indep", "cn", "ap", "so"),
    nq   = c(params$n_mean_qtl, params$n_shared_qtl, params$n_indep_qtl, params$n_cn_qtl, params$n_ap_qtl, n_so),
    va   = c(params$var_a_mean, params$var_a_shared, params$var_a_indep, params$var_a_cn, params$var_a_ap,
             if (is.null(params$var_a_so)) 0 else params$var_a_so),
    e2_b = c(0, params$ds_coupling^2 * rho^2 * S2, tau^2, rho^2 * phi_cn^2,
             params$ap_coupling^2 * rho^2 * S2, if (identical(so_axis, "driver")) rho^2 * S2 else 1),
    e2_a = c(1, 1 + params$ds_coupling^2 * s_h^2 / m, 0, p_cn^2 + p_cn * (1 - p_cn) / m,
             1 + params$ap_coupling^2 * s_h^2 / m, 0),
    stringsAsFactors = FALSE)
  one <- function(S2) {
    cls <- classes(S2)
    D_b <- sum(cls$nq * cls$va * cls$e2_b) + params$sigma_res^2 / (nr * (m - 1))
    D_a <- sum(cls$nq * cls$va * cls$e2_a) + params$sigma_res^2 / (nr * m)
    cls$ncp2_b <- n * cls$va * cls$e2_b / pmax(D_b - cls$va * cls$e2_b, 1e-12)
    cls$ncp2_a <- n * cls$va * cls$e2_a / pmax(D_a - cls$va * cls$e2_a, 1e-12)
    cls$power_slope <- mapply(epw, cls$va, cls$e2_b, MoreArgs = list(D = D_b))
    cls$power_int   <- mapply(epw, cls$va, cls$e2_a, MoreArgs = list(D = D_a))
    cls$D_slope <- D_b; cls$D_int <- D_a
    cls
  }
  out <- one(s_h^2)                       # quantities at E[S_h^2] = sigma_h^2
  if (average_env) {
    # S_h^2 ~ sigma_h^2 * chi-square(m - 1) / (m - 1); average power over it
    q  <- (seq_len(n_env_quad) - 0.5) / n_env_quad
    S2 <- s_h^2 * qchisq(q, m - 1) / (m - 1)
    pw_b <- rowMeans(vapply(S2, function(v) one(v)$power_slope, numeric(nrow(out))))
    pw_a <- rowMeans(vapply(S2, function(v) one(v)$power_int, numeric(nrow(out))))
    out$power_slope <- pw_b; out$power_int <- pw_a
  }
  out$tau <- tau; out$rho <- rho
  out[out$nq > 0, , drop = FALSE]
}

.SO_COLS    <- c("n_so_qtl", "so_axis", "so_center", "power_so_in_slope", "power_so_in_int",
                 "power_so_slope_only", "fdr_power_so", "win_power_so",
                 "so_slope_per_unit", "so_int_per_unit", "so_gmean_per_unit")
.CORE_COLS  <- c("core_sim", "core_maf", "fixed_effects")
.CLASSES    <- c("mean", "shared", "indep", "cn", "ap", "so")
.PRED_COLS  <- c("pred_mean_x", "pred_D_slope", "obs_var_slope", "pred_D_int", "obs_var_int",
                 paste0("pred_power_", .CLASSES, "_slope"), paste0("pred_power_", .CLASSES, "_int"),
                 paste0("pred_ncp2_", .CLASSES, "_slope"),
                 paste0("obs_power_", .CLASSES, "_int"))
.PVE_COLS   <- as.vector(outer(c("pve_env_", "pve_gmean_", "pve_slope_", "pve_int_"), .CLASSES, paste0))
.AUDIT_COLS <- c("cor_true_est_slope_bx", "H2_slope_bx", "H2_ratio_bx", "cor_true_est_slope_x",
                 "H2_slope_x", "n_sig_ratio_int_slope")
.v10_optional_cols <- c(.SO_COLS, .CORE_COLS, .PRED_COLS, .PVE_COLS, .AUDIT_COLS)

v10_optional_metrics <- function(gi, gs, rn, geno, pheno, idx, params, m_base) {
  out <- setNames(as.list(rep(NA_real_, length(.v10_optional_cols))), .v10_optional_cols)
  qi <- geno$qtl_idx; th <- params$alpha_bonf / params$n_markers
  si <- gi$marker[gi$p_value < th]; ss <- gs$marker[gs$p_value < th]
  pw <- function(q, s) if (length(q) > 0) mean(q %in% s) else NA_real_
  x <- idx$x; xc <- x - mean(x)

  if (length(qi$so) > 0) {
    v <- pheno$so_vec
    out$n_so_qtl <- length(qi$so); out$so_axis <- as.character(params$so_axis)
    out$so_center <- as.character(params$so_center)
    out$power_so_in_slope <- pw(qi$so, ss); out$power_so_in_int <- pw(qi$so, si)
    out$power_so_slope_only <- mean(qi$so %in% ss & !(qi$so %in% si))
    out$fdr_power_so <- pw(qi$so, gs$marker[gs$fdr < 0.05])
    out$win_power_so <- mean(vapply(qi$so, function(q) any(abs(ss - q) <= params$ld_window), logical(1)))
    b_v <- sum(xc * v) / sum(xc^2)
    out$so_slope_per_unit <- b_v                       # slope on x per unit g_so
    out$so_int_per_unit   <- mean(v) - b_v * mean(x)   # intercept at x = 0 per unit g_so
    out$so_gmean_per_unit <- mean(v)                   # environment-averaged mean per unit g_so
  }
  if (isTRUE(params$core_sim)) {
    out$core_sim <- TRUE; out$core_maf <- params$core_maf
    out$fixed_effects <- isTRUE(params$fixed_effects)
  } else if (isTRUE(params$fixed_effects)) out$fixed_effects <- TRUE

  need_ap <- isTRUE(params$predict_power) || isTRUE(params$pve_cols)
  if (need_ap) {
    eff <- locus_effect_matrix(geno, pheno, params)
    ap  <- analytic_power_conditional(geno, pheno, idx, params, eff)
    ap$A_full <- eff$A
  }
  if (isTRUE(params$predict_power)) {
    vpop <- function(u) mean((u - mean(u))^2)
    out$pred_mean_x <- ap$mean_x
    out$pred_D_slope <- ap$D_slope; out$obs_var_slope <- vpop(rn$slope)
    out$pred_D_int <- ap$D_int;     out$obs_var_int <- vpop(rn$intercept)
    for (cl in .CLASSES) {
      k <- ap$cls == cl
      if (!any(k)) next
      out[[paste0("pred_power_", cl, "_slope")]] <- mean(ap$power_slope[k])
      out[[paste0("pred_power_", cl, "_int")]]   <- mean(ap$power_int[k])
      out[[paste0("pred_ncp2_", cl, "_slope")]]  <- mean(ap$ncp2_slope[k])
      out[[paste0("obs_power_", cl, "_int")]]    <- pw(ap$pos[k], si)
    }
  }
  if (isTRUE(params$pve_cols)) {
    pv <- locus_pve(ap, pheno, rn)
    for (cl in unique(pv$cls)) {
      k <- pv$cls == cl
      out[[paste0("pve_env_", cl)]]   <- mean(pv$pve_env[k])
      out[[paste0("pve_gmean_", cl)]] <- mean(pv$pve_gmean[k])
      out[[paste0("pve_slope_", cl)]] <- mean(pv$pve_slope[k])
      out[[paste0("pve_int_", cl)]]   <- mean(pv$pve_int[k])
    }
  }
  if (isTRUE(params$audit_cols)) {
    # As in the Drive audit: slope scored against rho*true_slope_h + tau*true_slope_z.
    true_bx <- idx$actual_rho * pheno$true_slope_h + idx$tau * pheno$true_slope_z
    out$cor_true_est_slope_bx <- tryCatch(cor(true_bx, rn$slope, use = "complete.obs"), error = function(e) NA_real_)
    out$H2_slope_bx <- out$cor_true_est_slope_bx^2
    out$H2_ratio_bx <- if (!is.na(out$H2_slope_bx) && m_base$H2_int > 0) out$H2_slope_bx / m_base$H2_int else NA_real_
    # Exact: the true genotypic slope on x, sum_j c_j G_ij, from the generating values.
    cvec <- xc / sum(xc^2); m_e <- length(x)
    G <- outer(pheno$true_intercept, rep(1, m_e)) + outer(pheno$true_slope_h, pheno$h) +
      outer(pheno$true_slope_z, pheno$z) + outer(pheno$true_cn, pheno$cn_indicator) +
      outer(pheno$true_quad, pheno$h^2)
    if (!is.null(pheno$g_so)) G <- G + outer(pheno$g_so, pheno$so_vec)
    true_x <- as.vector(G %*% cvec)
    out$cor_true_est_slope_x <- tryCatch(cor(true_x, rn$slope, use = "complete.obs"), error = function(e) NA_real_)
    out$H2_slope_x <- out$cor_true_est_slope_x^2
    out$n_sig_ratio_int_slope <- if (m_base$n_sig_slope > 0) m_base$n_sig_int / m_base$n_sig_slope else NA_real_
  }
  out
}

# ============================================================================
# SECTION 8: VECTORIZED ROTATION TEST (Massive Speedup)
# ============================================================================

rotation_test <- function(geno, rn, rho_obs, n_perm = 200, alpha_bonf = 0.05) {
  X <- geno$X; n <- nrow(X); nm <- ncol(X); th <- alpha_bonf / nm
  gi_obs <- fast_gwas(X, rn$intercept)
  gs_obs <- fast_gwas(X, rn$slope)
  si_obs <- gi_obs$marker[gi_obs$p_value < th]
  ss_obs <- gs_obs$marker[gs_obs$p_value < th]
  unique_obs <- length(setdiff(ss_obs, si_obs))
  
  noise_sd <- sqrt(max(1 - rho_obs^2, 1e-10))
  
  Y_null <- matrix(rnorm(n * n_perm, 0, noise_sd), n, n_perm)
  Y_null <- sweep(Y_null, 1, scale(rn$intercept, scale = FALSE)[, 1] * rho_obs, "+")
  Y_null_c <- scale(Y_null, center = TRUE, scale = FALSE)
  
  ssy <- colSums(Y_null_c^2)
  ssx <- colSums(X^2)
  mono <- ssx < .Machine$double.eps; ssx[mono] <- 1
  
  XY <- crossprod(X, Y_null_c)
  beta <- XY / ssx
  ssr <- sweep(-(beta * XY), 2, ssy, "+")
  ssr[ssr < .Machine$double.eps] <- .Machine$double.eps
  df_res <- n - 2
  se <- sqrt(sweep(ssr / df_res, 1, ssx, "/"))
  tt <- beta / se
  pv <- 2 * pt(abs(tt), df = df_res, lower.tail = FALSE)
  pv[mono, ] <- 1
  
  if (length(si_obs) > 0) {
    null_unique <- colSums(pv[-si_obs, , drop = FALSE] < th)
  } else {
    null_unique <- colSums(pv < th)
  }
  
  list(unique_obs = unique_obs, null_mean = mean(null_unique), null_sd = sd(null_unique),
       null_q95 = quantile(null_unique, 0.95), null_q99 = quantile(null_unique, 0.99),
       pvalue = (sum(null_unique >= unique_obs) + 1) / (n_perm + 1), rho = rho_obs)
}

# ============================================================================
# SECTION 9: MULTI-INDEX GWAS
# ============================================================================

multi_index_gwas <- function(geno, pheno, idx_obj, alpha_bonf = 0.05) {
  X  <- geno$X; nm <- ncol(X); th <- alpha_bonf / nm; qi <- geno$qtl_idx
  
  rn_x <- estimate_reaction_norms(pheno$Y, idx_obj$x)
  gs_x <- fast_gwas(X, rn_x$slope)
  
  rn_z <- estimate_reaction_norms(pheno$Y, pheno$z)
  gs_z <- fast_gwas(X, rn_z$slope)
  
  ss_x <- gs_x$marker[gs_x$p_value < th]; ss_z <- gs_z$marker[gs_z$p_value < th]
  pw <- function(q, s) if (length(q) > 0) mean(q %in% s) else NA_real_
  
  list(
    n_sig_x = length(ss_x), n_sig_z = length(ss_z),
    overlap_xz = length(intersect(ss_x, ss_z)),
    unique_x = length(setdiff(ss_x, ss_z)), unique_z = length(setdiff(ss_z, ss_x)),
    power_shared_x = pw(qi$shared, ss_x), power_shared_z = pw(qi$shared, ss_z),
    power_indep_x  = pw(qi$indep,  ss_x), power_indep_z  = pw(qi$indep,  ss_z),
    power_cn_x = pw(qi$cn, ss_x), power_cn_z = pw(qi$cn, ss_z),
    power_ap_x = pw(qi$ap, ss_x), power_ap_z = pw(qi$ap, ss_z),
    cor_logp_xz = cor(-log10(gs_x$p_value), -log10(gs_z$p_value), method = "spearman")
  )
}

# ============================================================================
# SECTION 10: SINGLE REPLICATE RUNNER
# ============================================================================

run_replicate <- function(params, seed) {
  set.seed(seed)
  geno <- if (isTRUE(params$core_sim)) {
    simulate_genotypes_core(params$n_geno, params$n_markers, params$core_maf,
                            params$n_mean_qtl, params$n_shared_qtl, params$n_indep_qtl,
                            params$n_cn_qtl, params$n_ap_qtl)
  } else {
    simulate_genotypes(params$n_geno, params$n_markers, params$maf_min,
                       params$n_mean_qtl, params$n_shared_qtl, params$n_indep_qtl,
                       params$n_cn_qtl, params$n_ap_qtl, params$pop_structure)
  }
  pheno <- generate_phenotypes(geno, params)
  idx <- construct_index(pheno, params$rho, params$kappa)
  # v10: no-intercept slope QTL (no-op, and no RNG draw, when n_so_qtl = 0)
  so  <- add_slope_only_qtl(geno, pheno, idx, params, seed)
  geno <- so$geno; pheno <- so$pheno
  rn <- estimate_reaction_norms(pheno$Y, idx$x)
  
  na_slp <- is.na(rn$slope); na_int <- is.na(rn$intercept)
  if (any(na_slp)) rn$slope[na_slp] <- mean(rn$slope, na.rm = TRUE)
  if (any(na_int)) rn$intercept[na_int] <- mean(rn$intercept, na.rm = TRUE)
  if (any(is.na(rn$mean_perf))) rn$mean_perf[is.na(rn$mean_perf)] <- mean(rn$mean_perf, na.rm = TRUE)
  
  if (isTRUE(params$gwas_method == "pca")) {
    sub_idx <- sort(sample(ncol(geno$X), min(ncol(geno$X), 500)))
    sv <- svd(geno$X[, sub_idx, drop = FALSE], nu = params$n_pc, nv = 0)
    PCs <- sv$u[, seq_len(min(params$n_pc, ncol(sv$u))), drop = FALSE]
    gwas_fn <- function(X, y) fast_gwas_pca(X, y, n_pc = params$n_pc, PCs = PCs)
  } else if (isTRUE(params$gwas_method == "lmm")) {
    K <- tcrossprod(geno$X) / ncol(geno$X)
    gwas_fn <- function(X, y) fast_gwas_lmm(X, y, K = K)
  } else { 
    gwas_fn <- fast_gwas 
  }
  
  gi <- gwas_fn(geno$X, rn$intercept)
  gs <- gwas_fn(geno$X, rn$slope)
  
  m <- compute_metrics(gi, gs, rn, geno$qtl_idx, params$n_markers, pheno$Y, params$alpha_bonf, params$ld_window)
  
  m$rho_target <- round(params$rho, 4); m$rho_actual <- idx$actual_rho
  m$n_geno <- params$n_geno; m$n_env <- params$n_env
  m$n_indep_qtl <- params$n_indep_qtl; m$ds_coupling <- params$ds_coupling
  m$pop_structure <- params$pop_structure; m$seed <- seed
  
  m$cor_true_est_int <- cor(pheno$true_intercept, rn$intercept, use = "complete.obs")
  m$cor_true_est_slope <- tryCatch(cor(pheno$true_slope_h, rn$slope, use = "complete.obs"), error = function(e) NA)
  m$H2_int <- m$cor_true_est_int^2; m$H2_slope <- if (!is.na(m$cor_true_est_slope)) m$cor_true_est_slope^2 else NA_real_
  
  grand <- mean(pheno$Y, na.rm = TRUE); g_mns <- rowMeans(pheno$Y, na.rm = TRUE); e_mns <- colMeans(pheno$Y, na.rm = TRUE)
  MS_G <- (params$n_env * sum((g_mns - grand)^2)) / (params$n_geno - 1)
  MS_GxE <- (sum((pheno$Y - grand)^2, na.rm = TRUE) - (params$n_env * sum((g_mns - grand)^2)) - (params$n_geno * sum((e_mns - grand)^2))) / ((params$n_geno - 1) * (params$n_env - 1))
  m$H2_entry <- max((MS_G - MS_GxE) / MS_G, 0)
  m$H2_ratio <- if (!is.na(m$H2_slope) && m$H2_int > 0) m$H2_slope / m$H2_int else NA_real_
  
  qg <- compute_qg_metrics(pheno, idx, rn, geno, params)
  for (nm_q in names(qg)) m[[nm_q]] <- qg[[nm_q]]
  
  m$R_b_gwas_std <- if (!is.na(qg$R_b_ind_std)) qg$R_b_ind_std * sqrt(if (is.na(m$unique_prop)) 0 else m$unique_prop) else NA_real_
  
  # Genomic prediction CV (controlled by run_gp flag)
  if (isTRUE(params$run_gp)) {
    gp <- tryCatch(genomic_prediction_cv(geno$X, rn$intercept, rn$slope, 5), 
                   error = function(e) list(gp_acc_int = NA, gp_acc_slope = NA, gp_acc_int2 = NA, gp_acc_slope2 = NA))
  } else {
    gp <- list(gp_acc_int = NA, gp_acc_slope = NA, gp_acc_int2 = NA, gp_acc_slope2 = NA)
  }
  m$gp_acc_int <- gp$gp_acc_int; m$gp_acc_slope <- gp$gp_acc_slope
  m$gp_ratio <- if (!is.na(gp$gp_acc_slope2) && !is.na(gp$gp_acc_int2) && gp$gp_acc_int2 > 0.01) gp$gp_acc_slope2 / gp$gp_acc_int2 else NA_real_
  
  ammi <- tryCatch(estimate_ammi(pheno$Y), error = function(e) NULL)
  if (!is.null(ammi)) {
    m$cor_fwr_ammi <- cor(rn$slope, ammi$ipca1, use = "complete.obs")
    m$ammi_pct_ipca1 <- ammi$pct_ipca1[1]
  } else { m$cor_fwr_ammi <- m$ammi_pct_ipca1 <- NA_real_ }
  
  # v10 optional columns, appended after every v9 column. All NA when the v10
  # options are off; run_experiment() drops all-NA optional columns, so default
  # experiment tables keep the v9.0.0 column set.
  m <- c(m, v10_optional_metrics(gi, gs, rn, geno, pheno, idx, params, m))
  
  as.data.frame(m, stringsAsFactors = FALSE)
}

# ============================================================================
# SECTION 11: EXPERIMENT SUITE & PARALLELIZATION (ALL 25 EXPERIMENTS)
# ============================================================================

# --- List of function names to export to workers (no accumulated data) ---
.worker_fns <- c(
  "run_replicate", "simulate_genotypes", "generate_phenotypes",
  "construct_index", "estimate_reaction_norms",
  "fast_gwas", "fast_gwas_pca", "fast_gwas_lmm",
  "compute_metrics", "compute_qg_metrics",
  "genomic_prediction_cv", "estimate_ammi",
  "default_params",
  # v10
  "simulate_genotypes_core", "with_substream", "so_axis_vector", "add_slope_only_qtl",
  ".SO_SEED_OFFSET", "locus_effect_matrix", "analytic_power_conditional", "locus_pve",
  "v10_optional_metrics", ".v10_optional_cols", ".CLASSES", "closed_form_power"
)

run_experiment <- function(design, base_params, n_reps, seed_offset, label = "Experiment", n_cores = N_CORES) {
  n_conditions <- nrow(design); total <- n_conditions * n_reps
  cat(sprintf("  %s: %d conditions x %d reps = %d runs [%d workers]\n", label, n_conditions, n_reps, total, n_cores))
  
  # --- Build batched tasks (10 reps per task to reduce socket overhead) ---
  BATCH_SIZE <- 10L
  tasks <- list(); k <- 0L
  for (i in seq_len(n_conditions)) {
    p <- base_params
    for (nm in names(design)) p[[nm]] <- design[[nm]][i]
    # v10: design factors that run_replicate() does not already record (for example
    # kappa, env_range, gwas_method, sigma_res) are attached to each output row.
    # rho is recorded as rho_target and is not duplicated.
    extra <- lapply(setdiff(names(design), "rho"), function(nm) {
      v <- design[[nm]][i]; if (is.factor(v)) as.character(v) else v })
    names(extra) <- setdiff(names(design), "rho")
    seeds_all <- seed_offset + ((i - 1L) * n_reps) + seq_len(n_reps)
    batches <- split(seeds_all, ceiling(seq_along(seeds_all) / BATCH_SIZE))
    for (b in batches) {
      k <- k + 1L
      tasks[[k]] <- list(params = p, seeds = b, extra = extra)
    }
  }
  run_task <- function(task) {
    do.call(rbind, lapply(task$seeds, function(s) {
      r <- run_replicate(task$params, s)
      for (nm in names(task$extra)) {
        # attach when absent; also fill an optional v10 column left NA (e.g. n_so_qtl = 0)
        if (!(nm %in% names(r)) || (nm %in% .v10_optional_cols && is.na(r[[nm]][1]))) r[[nm]] <- task$extra[[nm]]
      }
      r
    }))
  }
  environment(run_task) <- globalenv()   # ship only the function, not this frame
  n_tasks <- length(tasks)
  
  t0 <- proc.time()
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Selective export: only functions and default_params, not accumulated results
    parallel::clusterExport(cl, .worker_fns, envir = globalenv())
    
    # Granular chunking for progress feedback
    n_chunks <- min(max(20L, ceiling(n_tasks / 50)), n_tasks)
    # cut() fails for a single interval (one task, as in Experiments 22-23 in
    # SMOKE_TEST mode), so a single chunk is used directly.
    chunk_indices <- if (n_chunks > 1L) {
      split(seq_len(n_tasks), cut(seq_len(n_tasks), n_chunks, labels = FALSE))
    } else list(seq_len(n_tasks))
    
    results_batched <- vector("list", n_tasks)
    for (c_idx in seq_len(n_chunks)) {
      cat(sprintf("    Chunk %d/%d (%d batches)...\n", c_idx, n_chunks, length(chunk_indices[[c_idx]])))
      chunk_res <- parallel::parLapply(cl, tasks[chunk_indices[[c_idx]]], run_task)
      results_batched[chunk_indices[[c_idx]]] <- chunk_res
    }
    
    # Flatten batched results
    results <- do.call(rbind, results_batched)
    
  } else {
    results <- do.call(rbind, lapply(tasks, run_task))
  }
  
  # v10: drop optional columns that are NA in every row (v10 options off), so
  # default experiments keep exactly the v9.0.0 column set.
  opt <- intersect(names(results), .v10_optional_cols)
  drop <- opt[vapply(opt, function(cl) all(is.na(results[[cl]])), logical(1))]
  if (length(drop)) results <- results[, setdiff(names(results), drop), drop = FALSE]
  
  cat(sprintf("    Complete. %.1fs\n", (proc.time() - t0)[3]))
  results
}

# 1-11 Core & Architectures
run_rho_sweep <- function(n_reps = 1000) run_experiment(data.frame(rho = rho_grid_fine), default_params, n_reps, 20260221L, "Rho sweep")
run_null_model <- function(n_reps = 300) { p <- default_params; p$n_indep_qtl <- 0; p$n_cn_qtl <- 0; p$n_ap_qtl <- 0; run_experiment(data.frame(rho = rho_grid_fine), p, n_reps, 40000000L, "Null model") }
run_kappa_comparison <- function(n_reps = 200) run_experiment(expand.grid(rho = rho_grid_coarse, kappa = c(0, 0.25, 0.5, 0.75, 1)), default_params, n_reps, 30000000L, "Kappa comparison")
run_ds_coupling_sweep <- function(n_reps = 200) run_experiment(expand.grid(rho = rho_grid_coarse, ds_coupling = c(0, 0.3, 0.6, 0.9, 1)), default_params, n_reps, 50000000L, "DS coupling sweep")
run_env_range_sweep <- function(n_reps = 200) run_experiment(expand.grid(rho = rho_grid_coarse, env_range = c("wide", "narrow")), default_params, n_reps, 55000000L, "Env range sweep")
run_nenv_sweep <- function(n_reps = 200) run_experiment(expand.grid(rho = rho_grid_coarse, n_env = c(3, 4, 5, 6, 8, 12, 16, 20)), default_params, n_reps, 60000000L, "N-env sweep")
run_pop_structure_contrast <- function(n_reps = 150) run_experiment(expand.grid(rho = rho_grid_coarse, pop_structure = c("diverse", "biparental", "nam")), default_params, n_reps, 65000000L, "Pop structure")

run_architecture_comparison <- function(n_reps = 200) {
  configs <- list(
    ds_only  = list(n_shared_qtl = 10, n_indep_qtl = 0, n_cn_qtl = 0, n_ap_qtl = 0),
    ds_indep = list(n_shared_qtl = 10, n_indep_qtl = 5, n_cn_qtl = 0, n_ap_qtl = 0),
    ds_cn    = list(n_shared_qtl = 10, n_indep_qtl = 0, n_cn_qtl = 5, n_ap_qtl = 0),
    ds_ap    = list(n_shared_qtl = 10, n_indep_qtl = 0, n_cn_qtl = 0, n_ap_qtl = 5),
    full     = list(n_shared_qtl = 10, n_indep_qtl = 5, n_cn_qtl = 3, n_ap_qtl = 2)
  )
  all_results <- vector("list", length(configs))
  for (ci in seq_along(configs)) {
    p <- default_params; for (nm in names(configs[[ci]])) p[[nm]] <- configs[[ci]][[nm]]
    res <- run_experiment(data.frame(rho = rho_grid_coarse), p, n_reps, 70000000L + ci * 100000L, paste("Architecture:", names(configs)[ci]))
    res$architecture <- names(configs)[ci]; all_results[[ci]] <- res
  }
  do.call(rbind, all_results)
}

run_sample_size_sweep <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, n_geno = c(50, 100, 200, 500, 1000, 2000, 5000)), default_params, n_reps, 75000000L, "Genotype size")
run_nonlinear_test <- function(n_reps = 200) { p <- default_params; p$var_a_quad <- 0.03; run_experiment(expand.grid(rho = rho_grid_coarse, nonlinear = c(FALSE, TRUE)), p, n_reps, 80000000L, "Nonlinear test") }
run_effect_size_sensitivity <- function(n_reps = 200) run_experiment(expand.grid(rho = c(0.00, 0.50, 0.85, 0.95, 0.996), var_a_indep = c(0.01, 0.02, 0.04, 0.08)), default_params, n_reps, 85000000L, "Effect size sensitivity")

# 12 Rotation test
run_rotation_test_validation <- function(n_reps = 50, n_perm = 200) {
  cat(sprintf("  Rotation test: %d rho levels x %d reps x %d perms\n", length(rho_grid_coarse), n_reps, n_perm))
  tasks <- expand.grid(rho = rho_grid_coarse, rep = 1:n_reps)
  cl <- parallel::makeCluster(N_CORES); on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c(.worker_fns, "rotation_test"), envir = globalenv())
  
  # Batch rotation tasks too
  task_list <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / 10L))
  res_list <- vector("list", length(task_list))
  for (ci in seq_along(task_list)) {
    if (ci %% 5 == 1) cat(sprintf("    Rotation chunk %d/%d...\n", ci, length(task_list)))
    res_list[[ci]] <- parallel::parLapply(cl, task_list[[ci]], function(i) {
      rv <- tasks$rho[i]
      set.seed(90000000L + i)
      p <- default_params; p$rho <- rv
      g <- simulate_genotypes(p$n_geno, p$n_markers, p$maf_min, p$n_mean_qtl, p$n_shared_qtl, p$n_indep_qtl, p$n_cn_qtl, p$n_ap_qtl, p$pop_structure)
      ph <- generate_phenotypes(g, p); ix <- construct_index(ph, rv, 1.0)
      rn <- estimate_reaction_norms(ph$Y, ix$x)
      rt <- rotation_test(g, rn, ix$actual_rho, n_perm, p$alpha_bonf)
      data.frame(rho_target = rv, rho_actual = ix$actual_rho, unique_obs = rt$unique_obs, null_mean = rt$null_mean, pvalue = rt$pvalue)
    })
  }
  do.call(rbind, unlist(res_list, recursive = FALSE))
}

# 13 Multi-index
run_multi_index_comparison <- function(n_reps = 150) {
  cat(sprintf("  Multi-index: %d rho levels x %d reps\n", length(rho_grid_coarse), n_reps))
  tasks <- expand.grid(rho = rho_grid_coarse, rep = 1:n_reps)
  cl <- parallel::makeCluster(N_CORES); on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c(.worker_fns, "multi_index_gwas"), envir = globalenv())
  
  task_list <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / 10L))
  res_list <- vector("list", length(task_list))
  for (ci in seq_along(task_list)) {
    if (ci %% 5 == 1) cat(sprintf("    Multi-index chunk %d/%d...\n", ci, length(task_list)))
    res_list[[ci]] <- parallel::parLapply(cl, task_list[[ci]], function(i) {
      rv <- tasks$rho[i]
      set.seed(95000000L + i)
      p <- default_params; p$rho <- rv
      g <- simulate_genotypes(p$n_geno, p$n_markers, p$maf_min, p$n_mean_qtl, p$n_shared_qtl, p$n_indep_qtl, p$n_cn_qtl, p$n_ap_qtl, p$pop_structure)
      ph <- generate_phenotypes(g, p); ix <- construct_index(ph, rv, 1.0)
      mi <- multi_index_gwas(g, ph, ix, p$alpha_bonf)
      mi$rho_target <- round(rv, 4); mi$rho_actual <- ix$actual_rho; mi$seed <- 95000000L + i
      as.data.frame(mi, stringsAsFactors = FALSE)
    })
  }
  do.call(rbind, unlist(res_list, recursive = FALSE))
}

# 14-25 Extensions, LMM, & Allele Frequencies
run_marker_density_sweep <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, n_markers = c(500, 1000, 2000, 5000, 10000)), default_params, n_reps, 95000000L, "Marker density")
run_heritability_sweep <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, sigma_res = c(0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 8.0)), default_params, n_reps, 100000000L, "Heritability sweep")

run_qtl_count_sweep <- function(n_reps = 100) {
  qtl_configs <- list(
    list(nm = 3, ns = 3, ni = 2, nc = 1, na = 1), 
    list(nm = 10, ns = 10, ni = 5, nc = 3, na = 2), 
    list(nm = 17, ns = 17, ni = 8, nc = 5, na = 3), 
    list(nm = 33, ns = 33, ni = 17, nc = 10, na = 7),
    list(nm = 165, ns = 165, ni = 85, nc = 50, na = 35),
    list(nm = 330, ns = 330, ni = 170, nc = 100, na = 70)
  )
  all_results <- vector("list", length(qtl_configs))
  for (ci in seq_along(qtl_configs)) {
    cfg <- qtl_configs[[ci]]
    n_total <- cfg$nm + cfg$ns + cfg$ni + cfg$nc + cfg$na
    base_var <- 0.04 * 30 / n_total
    p <- default_params
    p$n_mean_qtl <- cfg$nm; p$n_shared_qtl <- cfg$ns; p$n_indep_qtl <- cfg$ni; p$n_cn_qtl <- cfg$nc; p$n_ap_qtl <- cfg$na
    p$var_a_mean <- base_var; p$var_a_shared <- base_var; p$var_a_indep <- base_var; p$var_a_cn <- base_var; p$var_a_ap <- base_var
    res <- run_experiment(data.frame(rho = rho_grid_sensitivity), p, n_reps, 105000000L + ci * 10000L, sprintf("QTL count: %d total", n_total))
    res$n_total_qtl <- n_total; all_results[[ci]] <- res
  }
  do.call(rbind, all_results)
}

run_heterogeneous_residuals <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, hetero_res = c(FALSE, TRUE)), default_params, n_reps, 110000000L, "Heterogeneous residuals")
run_pca_gwas_robustness <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, gwas_method = c("simple", "pca", "lmm"), pop_structure = c("diverse", "nam")), default_params, n_reps, 115000000L, "PCA/LMM Robustness")

run_effect_distribution <- function(n_reps = 100) {
  res_eq <- run_experiment(data.frame(rho = rho_grid_sensitivity), default_params, n_reps, 120000000L, "Effect dist: equal")
  res_eq$effect_dist <- "equal"
  p_exp <- default_params; p_exp$var_a_shared <- 0.10; p_exp$var_a_indep <- 0.01; p_exp$var_a_cn <- 0.01; p_exp$var_a_ap <- 0.02
  res_exp <- run_experiment(data.frame(rho = rho_grid_sensitivity), p_exp, n_reps, 121000000L, "Effect dist: unequal")
  res_exp$effect_dist <- "unequal"; rbind(res_eq, res_exp)
}

run_missing_data_sweep <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, pct_missing = c(0, 0.05, 0.1, 0.2)), default_params, n_reps, 125000000L, "Missing data")
run_replication_sweep <- function(n_reps = 100) run_experiment(expand.grid(rho = rho_grid_sensitivity, n_reps = c(1L, 2L, 3L, 4L)), default_params, n_reps, 130000000L, "Replication sweep")

run_convergence_diagnostic <- function(rep_scale = 1) {
  cat("    Running convergence diagnostic...\n")
  rep_levels <- sort(unique(pmax(2L, as.integer(round(c(10, 25, 50, 100, 150, 200, 300) * rep_scale)))))
  rho_test <- c(0.50, 0.85, 0.996)
  all_results <- vector("list", length(rep_levels) * length(rho_test))
  k <- 0L
  for (rv in rho_test) {
    p <- default_params; p$rho <- rv
    res_full <- run_experiment(data.frame(rho = rv), p, max(rep_levels), 135000000L + round(rv * 1000), sprintf("Convergence: rho=%.3f", rv))
    for (nr in rep_levels) {
      k <- k + 1L; sub <- res_full[seq_len(nr), ]
      all_results[[k]] <- data.frame(rho_target = rv, n_reps_used = nr, mean_overlap = mean(sub$overlap_prop, na.rm = TRUE), se_overlap = sd(sub$overlap_prop, na.rm = TRUE) / sqrt(nr), mean_power_indep = mean(sub$power_indep_in_slope, na.rm = TRUE), se_power_indep = sd(sub$power_indep_in_slope, na.rm = TRUE) / sqrt(nr), mean_cor_si = mean(sub$cor_slope_int, na.rm = TRUE), se_cor_si = sd(sub$cor_slope_int, na.rm = TRUE) / sqrt(nr), mean_gp_slope = mean(sub$gp_acc_slope, na.rm = TRUE), se_gp_slope = sd(sub$gp_acc_slope, na.rm = TRUE) / sqrt(nr), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, all_results)
}

run_seed_sensitivity <- function(n_reps = 100) {
  seeds <- c(1000000L, 2000000L, 3000000L, 4000000L, 5000000L); rho_test <- c(0.50, 0.85, 0.996)
  all_results <- vector("list", length(seeds) * length(rho_test)); k <- 0L
  for (si in seq_along(seeds)) {
    for (rv in rho_test) {
      k <- k + 1L; p <- default_params; p$rho <- rv
      res <- run_experiment(data.frame(rho = rv), p, n_reps, seeds[si] + round(rv * 1000), sprintf("Seed %d, rho=%.3f", seeds[si], rv))
      all_results[[k]] <- data.frame(seed_group = si, rho_target = rv, mean_overlap = mean(res$overlap_prop, na.rm = TRUE), mean_power_indep = mean(res$power_indep_in_slope, na.rm = TRUE), mean_cor_si = mean(res$cor_slope_int, na.rm = TRUE), mean_gp_slope = mean(res$gp_acc_slope, na.rm = TRUE), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, all_results)
}

run_ld_window_sweep <- function(n_reps = 100) {
  all_results <- list()
  for (win in c(0L, 2L, 5L, 10L)) {
    p <- default_params; p$ld_window <- win
    res <- run_experiment(data.frame(rho = rho_grid_sensitivity), p, n_reps, 140000000L + win * 10000L, sprintf("LD window: %d", win))
    res$ld_window <- win; all_results <- c(all_results, list(res))
  }
  do.call(rbind, all_results)
}

run_maf_sweep <- function(n_reps = 100) {
  run_experiment(expand.grid(rho = rho_grid_sensitivity, maf_min = c(0.01, 0.05, 0.10, 0.20, 0.35)), default_params, n_reps, 150000000L, "Minor Allele Frequency (MAF) sweep")
}

# ============================================================================
# SECTION 11B: v10 EXPERIMENTS (26-29)
# ============================================================================

# 26 No-intercept slope QTL class. Default architecture plus 5 slope-only loci
# (var_a_so = 0.04, the per-locus variance of every other class). With kappa = 1,
# mean(x) = 0 and the two centerings coincide, so centering is varied at kappa = 0.5.
run_slope_only_qtl_sweep <- function(n_reps = 200) {
  p <- default_params; p$n_so_qtl <- 5L; p$predict_power <- TRUE
  design <- rbind(
    expand.grid(rho = rho_grid_coarse, kappa = 1.0, so_axis = c("index", "driver"),
                so_center = "intercept", stringsAsFactors = FALSE),
    expand.grid(rho = rho_grid_coarse, kappa = 0.5, so_axis = c("index", "driver"),
                so_center = c("intercept", "mean"), stringsAsFactors = FALSE))
  run_experiment(design, p, n_reps, 160000000L, "Slope-only QTL class")
}

# 27 Simplified core simulation: independent loci, fixed allele frequency 0.30,
# fixed effect magnitudes, no LD; default class counts plus 5 index-axis slope-only loci.
run_core_simulation <- function(n_reps = 200) {
  p <- default_params; p$core_sim <- TRUE; p$fixed_effects <- TRUE
  p$n_so_qtl <- 5L; p$predict_power <- TRUE
  run_experiment(data.frame(rho = rho_grid_coarse), p, n_reps, 170000000L, "Core simulation")
}

# 28 Within-environment variance (core simulation, analytic prediction on).
#   28a: residual SD x within-environment replicates x mean-axis G x E present/absent
#   28b: main-effect genetic variance x environmental range x number of environments
run_within_env_variance <- function(n_reps = 100) {
  p <- default_params; p$core_sim <- TRUE; p$fixed_effects <- TRUE
  p$n_so_qtl <- 5L; p$predict_power <- TRUE; p$run_gp <- FALSE
  scen <- data.frame(gxe_mean_axis = c("default", "none"), ds_coupling = c(0.8, 0),
                     n_ap_qtl = c(2L, 0L), n_cn_qtl = c(3L, 0L), stringsAsFactors = FALSE)
  a <- expand.grid(rho = rho_grid_sensitivity, sigma_res = c(0.5, 1, 2, 4), n_reps = c(1L, 4L),
                   s = seq_len(nrow(scen)))
  a <- cbind(a[, c("rho", "sigma_res", "n_reps")], scen[a$s, ], part = "28a",
             var_a_mean = 0.04, env_range = "wide", n_env = 12)
  b <- expand.grid(rho = rho_grid_sensitivity, var_a_mean = c(0.04, 0.4),
                   env_range = c("wide", "narrow"), n_env = c(6, 12, 24), stringsAsFactors = FALSE)
  b <- cbind(b, sigma_res = 1, n_reps = 1L, scen[1, ], part = "28b")
  design <- rbind(a[, names(b)], b)
  rownames(design) <- NULL
  run_experiment(design, p, n_reps, 180000000L, "Within-environment variance")
}

# 29 Per-locus variance explained by QTL class at default parameters, with and
# without 5 index-axis slope-only loci.
run_locus_pve <- function(n_reps = 200) {
  p <- default_params; p$pve_cols <- TRUE; p$run_gp <- FALSE
  design <- expand.grid(rho = rho_grid_sensitivity, n_so_qtl = c(0L, 5L))
  run_experiment(design, p, n_reps, 190000000L, "Locus variance explained")
}

summarize_v10 <- function(res, by_vars, patterns = c("^power_", "^pred_", "^obs_", "^pve_",
                                                     "^so_", "^fdr_power", "^win_power",
                                                     "^cor_slope_int$", "^overlap_prop$",
                                                     "^unique_prop$", "^rho_actual$")) {
  cols <- unique(unlist(lapply(patterns, function(p) grep(p, names(res), value = TRUE))))
  cols <- cols[vapply(cols, function(cl) is.numeric(res[[cl]]), logical(1))]
  grp <- interaction(res[, by_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(res, grp), function(d) {
    row <- d[1, by_vars, drop = FALSE]; row$n_reps_run <- nrow(d)
    for (cl in cols) row[[cl]] <- mean(d[[cl]], na.rm = TRUE)
    row
  }))
  rownames(out) <- NULL
  out[do.call(order, unname(as.list(out[, by_vars, drop = FALSE]))), , drop = FALSE]
}

# ============================================================================
# SECTION 12: SUMMARY FUNCTIONS
# ============================================================================

summarize_rho_sweep <- function(res) {
  # Safe mean helper: returns NA if column doesn't exist
  sm <- function(d, col) if (col %in% names(d)) mean(d[[col]], na.rm = TRUE) else NA_real_
  
  do.call(rbind, lapply(split(res, res$rho_target), function(d) {
    data.frame(
      rho_target = d$rho_target[1], n_reps = nrow(d),
      rho_actual = sm(d, "rho_actual"),
      cor_slope_int = sm(d, "cor_slope_int"), cor_si_sd = sd(d$cor_slope_int, na.rm = TRUE),
      overlap_prop = sm(d, "overlap_prop"), unique_prop = sm(d, "unique_prop"),
      power_indep = sm(d, "power_indep_in_slope"), power_shared = sm(d, "power_shared_in_slope"),
      power_cn = sm(d, "power_cn_in_slope"), power_ap = sm(d, "power_ap_in_slope"),
      cor_logp = sm(d, "cor_logp"),
      fpr_slope = sm(d, "fpr_slope"),
      median_rk_cor = sm(d, "median_rank_cor"), prop_low_rk = sm(d, "prop_low_rank_cor"),
      H2_int = sm(d, "H2_int"), H2_slope = sm(d, "H2_slope"),
      H2_entry = sm(d, "H2_entry"), H2_ratio = sm(d, "H2_ratio"),
      r_predicted = sm(d, "r_predicted"), r_true_genetic = sm(d, "r_true_genetic"),
      noise_inflation = sm(d, "noise_inflation"),
      V_b_cond_a = sm(d, "V_b_cond_a"), pct_V_b_unique = sm(d, "pct_V_b_unique"),
      mean_type_b = sm(d, "mean_type_b"), min_type_b = sm(d, "min_type_b"),
      r_g_via_lande = sm(d, "r_g_via_lande"), n_eff_grn = sm(d, "n_eff_grn"),
      h2_b_narrow = sm(d, "h2_b_narrow"),
      R_b_std = sm(d, "R_b_std"), R_b_ind_std = sm(d, "R_b_ind_std"),
      R_b_gwas_std = sm(d, "R_b_gwas_std"),
      n_eff_env = sm(d, "n_eff_env"), pct_pc1 = sm(d, "pct_pc1"),
      fdr_overlap = sm(d, "fdr_overlap_prop"), fdr_unique = sm(d, "fdr_unique_prop"),
      fdr_pw_indep = sm(d, "fdr_power_indep"), fdr_pw_shared = sm(d, "fdr_power_shared"),
      lambda_int = sm(d, "lambda_int"), lambda_slope = sm(d, "lambda_slope"),
      gp_acc_int = sm(d, "gp_acc_int"), gp_acc_slope = sm(d, "gp_acc_slope"),
      gp_ratio = sm(d, "gp_ratio"),
      cor_fwr_ammi = sm(d, "cor_fwr_ammi"), ammi_pct_ipca1 = sm(d, "ammi_pct_ipca1"),
      win_overlap = sm(d, "win_overlap_prop"), win_unique = sm(d, "win_unique_prop"),
      win_pw_indep = sm(d, "win_power_indep")
    )
  }))
}

summarize_by <- function(res, by_vars) {
  sm <- function(d, col) if (col %in% names(d)) mean(d[[col]], na.rm = TRUE) else NA_real_
  
  grp <- interaction(res[, by_vars, drop = FALSE], drop = TRUE)
  do.call(rbind, lapply(split(res, grp), function(d) {
    row <- d[1, by_vars, drop = FALSE]; row$n_reps <- nrow(d)
    row$cor_slope_int <- sm(d, "cor_slope_int"); row$overlap_prop <- sm(d, "overlap_prop")
    row$unique_prop <- sm(d, "unique_prop"); row$power_indep <- sm(d, "power_indep_in_slope")
    row$power_shared <- sm(d, "power_shared_in_slope"); row$power_cn <- sm(d, "power_cn_in_slope")
    row$power_ap <- sm(d, "power_ap_in_slope"); row$cor_logp <- sm(d, "cor_logp")
    row$median_rk_cor <- sm(d, "median_rank_cor"); row$H2_int <- sm(d, "H2_int")
    row$H2_slope <- sm(d, "H2_slope"); row$H2_entry <- sm(d, "H2_entry")
    row$H2_ratio <- sm(d, "H2_ratio"); row$r_predicted <- sm(d, "r_predicted")
    row$r_true_genetic <- sm(d, "r_true_genetic"); row$noise_inflation <- sm(d, "noise_inflation")
    row$V_b_cond_a <- sm(d, "V_b_cond_a"); row$pct_V_b_unique <- sm(d, "pct_V_b_unique")
    row$mean_type_b <- sm(d, "mean_type_b"); row$min_type_b <- sm(d, "min_type_b")
    row$n_eff_grn <- sm(d, "n_eff_grn"); row$h2_b_narrow <- sm(d, "h2_b_narrow")
    row$R_b_std <- sm(d, "R_b_std"); row$R_b_ind_std <- sm(d, "R_b_ind_std")
    row$R_b_gwas_std <- sm(d, "R_b_gwas_std"); row$n_eff_env <- sm(d, "n_eff_env")
    row$pct_pc1 <- sm(d, "pct_pc1"); row$fdr_overlap <- sm(d, "fdr_overlap_prop")
    row$fdr_unique <- sm(d, "fdr_unique_prop"); row$fdr_pw_indep <- sm(d, "fdr_power_indep")
    row$fdr_pw_shared <- sm(d, "fdr_power_shared"); row$lambda_int <- sm(d, "lambda_int")
    row$lambda_slope <- sm(d, "lambda_slope"); row$gp_acc_int <- sm(d, "gp_acc_int")
    row$gp_acc_slope <- sm(d, "gp_acc_slope"); row$gp_ratio <- sm(d, "gp_ratio")
    row$cor_fwr_ammi <- sm(d, "cor_fwr_ammi"); row$ammi_pct_ipca1 <- sm(d, "ammi_pct_ipca1")
    row$win_overlap <- sm(d, "win_overlap_prop"); row$win_unique <- sm(d, "win_unique_prop")
    row$win_pw_indep <- sm(d, "win_power_indep"); row
  }))
}

# --- Bootstrap Inference ---

bc_ci <- function(boot_means, obs_mean, alpha = 0.05) {
  z0 <- qnorm(mean(boot_means <= obs_mean, na.rm = TRUE))
  if (!is.finite(z0)) return(quantile(boot_means, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  za <- qnorm(alpha / 2); lo_p <- max(min(pnorm(2 * z0 + za), 1 - 1e-6), 1e-6); hi_p <- max(min(pnorm(2 * z0 - za), 1 - 1e-6), 1e-6)
  quantile(boot_means, c(lo_p, hi_p), na.rm = TRUE)
}

boot_vector <- function(x, n_boot = 2000, alpha = 0.05) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 5) return(list(mean = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, se = NA_real_))
  obs_mean <- mean(x); boot_means <- replicate(n_boot, mean(sample(x, n, replace = TRUE)))
  ci <- bc_ci(boot_means, obs_mean, alpha)
  list(mean = obs_mean, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]), se = sd(boot_means))
}

bootstrap_summary <- function(res, by_var = "rho_target", metrics = NULL, n_boot = 2000, alpha = 0.05) {
  if (is.null(metrics)) metrics <- c("cor_slope_int", "overlap_prop", "unique_prop", "power_indep_in_slope", "power_shared_in_slope", "power_cn_in_slope", "power_ap_in_slope", "cor_logp", "H2_int", "H2_slope", "H2_entry", "H2_ratio", "fpr_slope", "cor_true_est_int", "cor_true_est_slope", "r_predicted", "r_true_genetic", "noise_inflation", "V_b_cond_a", "pct_V_b_unique", "mean_type_b", "min_type_b", "r_g_via_lande", "n_eff_grn", "h2_b_narrow", "R_b_std", "R_b_ind_std", "R_b_gwas_std", "n_eff_env", "pct_pc1")
  metrics <- intersect(metrics, names(res)); levels <- sort(unique(res[[by_var]]))
  out <- vector("list", length(levels) * length(metrics)); k <- 0L
  for (lv in levels) {
    sub <- res[res[[by_var]] == lv, ]
    for (met in metrics) { k <- k + 1L; bs <- boot_vector(sub[[met]], n_boot, alpha); out[[k]] <- data.frame(level = lv, metric = met, n = sum(!is.na(sub[[met]])), mean = bs$mean, ci_lo = bs$ci_lo, ci_hi = bs$ci_hi, se = bs$se, stringsAsFactors = FALSE, row.names = NULL) }
  }
  do.call(rbind, out)
}

boot_contrast_pair <- function(x1, x2, n_boot = 2000, alpha = 0.05) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 5 || length(x2) < 5) return(list(diff = NA, ci_lo = NA, ci_hi = NA, p_value = NA, cohens_d = NA, d_ci_lo = NA, d_ci_hi = NA))
  obs_diff <- mean(x2) - mean(x1); pooled_sd <- sqrt((var(x1) * (length(x1) - 1) + var(x2) * (length(x2) - 1)) / (length(x1) + length(x2) - 2))
  obs_d <- if (pooled_sd > 1e-10) obs_diff / pooled_sd else NA_real_
  boot_mat <- replicate(n_boot, { b1 <- sample(x1, length(x1), replace = TRUE); b2 <- sample(x2, length(x2), replace = TRUE); bd <- mean(b2) - mean(b1); bsd <- sqrt((var(b1) * (length(b1) - 1) + var(b2) * (length(b2) - 1)) / (length(b1) + length(b2) - 2)); c(diff = bd, d = if (bsd > 1e-10) bd / bsd else NA) })
  boot_diffs <- boot_mat["diff", ]; boot_ds <- boot_mat["d", ]
  p_val <- min(if (obs_diff >= 0) 2 * mean(boot_diffs <= 0) else 2 * mean(boot_diffs >= 0), 1.0)
  ci_diff <- bc_ci(boot_diffs, obs_diff, alpha); d_ci <- quantile(boot_ds, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  list(diff = obs_diff, ci_lo = unname(ci_diff[1]), ci_hi = unname(ci_diff[2]), p_value = p_val, cohens_d = obs_d, d_ci_lo = unname(d_ci[1]), d_ci_hi = unname(d_ci[2]))
}

bootstrap_contrasts <- function(res, by_var = "rho_target", pairs = NULL, metrics = NULL, n_boot = 2000, alpha = 0.05) {
  if (is.null(metrics)) metrics <- c("cor_slope_int", "overlap_prop", "unique_prop", "power_indep_in_slope", "power_shared_in_slope", "power_cn_in_slope", "power_ap_in_slope", "H2_int", "H2_slope", "H2_ratio", "pct_V_b_unique", "R_b_std", "R_b_ind_std", "R_b_gwas_std", "n_eff_env", "mean_type_b", "min_type_b", "noise_inflation")
  metrics <- intersect(metrics, names(res))
  if (is.null(pairs)) pairs <- list(c(0.00, 0.50), c(0.50, 0.85), c(0.85, 0.95), c(0.95, 0.996), c(0.00, 0.996))
  actual_levels <- sort(unique(res[[by_var]])); snap <- function(v) actual_levels[which.min(abs(actual_levels - v))]
  out <- vector("list", length(pairs) * length(metrics)); k <- 0L
  for (pr in pairs) {
    lv1 <- snap(pr[1]); lv2 <- snap(pr[2]); g1 <- res[res[[by_var]] == lv1, ]; g2 <- res[res[[by_var]] == lv2, ]
    for (met in metrics) { k <- k + 1L; bc <- boot_contrast_pair(g1[[met]], g2[[met]], n_boot, alpha); out[[k]] <- data.frame(level_1 = pr[1], level_2 = pr[2], metric = met, mean_1 = mean(g1[[met]], na.rm = TRUE), mean_2 = mean(g2[[met]], na.rm = TRUE), diff = bc$diff, ci_lo = bc$ci_lo, ci_hi = bc$ci_hi, p_value = bc$p_value, cohens_d = bc$cohens_d, d_ci_lo = bc$d_ci_lo, d_ci_hi = bc$d_ci_hi, stringsAsFactors = FALSE, row.names = NULL) }
  }
  result <- do.call(rbind, out); result$p_adjusted <- p.adjust(result$p_value, method = "holm"); result$sig <- ifelse(result$p_adjusted < 0.001, "***", ifelse(result$p_adjusted < 0.01, "**", ifelse(result$p_adjusted < 0.05, "*", "ns"))); result
}

bootstrap_architecture_contrasts <- function(res_arch, n_boot = 2000, alpha = 0.05) {
  metrics <- intersect(c("cor_slope_int", "overlap_prop", "unique_prop", "power_indep_in_slope", "power_cn_in_slope", "power_ap_in_slope", "H2_slope", "H2_ratio", "pct_V_b_unique", "R_b_std", "R_b_gwas_std", "n_eff_env"), names(res_arch))
  rho_levels <- sort(unique(res_arch$rho_target)); out <- vector("list", length(rho_levels) * length(metrics)); k <- 0L
  for (rv in rho_levels) {
    full_sub <- res_arch[res_arch$rho_target == rv & res_arch$architecture == "full", ]; null_sub <- res_arch[res_arch$rho_target == rv & res_arch$architecture == "ds_only", ]
    for (met in metrics) { k <- k + 1L; bc <- boot_contrast_pair(null_sub[[met]], full_sub[[met]], n_boot, alpha); out[[k]] <- data.frame(rho = rv, contrast = "full_vs_ds_only", metric = met, mean_ds_only = mean(null_sub[[met]], na.rm = TRUE), mean_full = mean(full_sub[[met]], na.rm = TRUE), diff = bc$diff, ci_lo = bc$ci_lo, ci_hi = bc$ci_hi, p_value = bc$p_value, cohens_d = bc$cohens_d, stringsAsFactors = FALSE, row.names = NULL) }
  }
  result <- do.call(rbind, out); result$p_adjusted <- p.adjust(result$p_value, method = "holm"); result$sig <- ifelse(result$p_adjusted < 0.001, "***", ifelse(result$p_adjusted < 0.01, "**", ifelse(result$p_adjusted < 0.05, "*", "ns"))); result
}

# ============================================================================
# SECTION 13: PUBLICATION FIGURES
# ============================================================================

pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#000000")
ag <- function(df, by, val) { s <- split(df[[val]], df[[by]]); data.frame(x = as.numeric(names(s)), mn = sapply(s, mean, na.rm = TRUE), lo = sapply(s, quantile, 0.025, na.rm = TRUE), hi = sapply(s, quantile, 0.975, na.rm = TRUE), row.names = NULL) }

plot_results <- function(res_rho, res_null = NULL, res_arch = NULL, res_ds = NULL, res_nenv = NULL, res_pop = NULL, res_rotation = NULL, res_multi = NULL, outdir = ".") {
  pdf(file.path(outdir, "rho_sweep_cor.pdf"), width = 5.5, height = 4.5); par(mar = c(4.8, 5, 1.5, 1), family = "serif", cex.lab = 1.3, cex.axis = 1.1)
  a1 <- ag(res_rho, "rho_target", "cor_slope_int"); plot(a1$x, a1$mn, type = "n", xlim = c(-0.03, 1.01), ylim = c(-0.1, 1.05), xlab = expression("Index-mean correlation (" * rho * ")"), ylab = expression("cor(slope, intercept)"))
  polygon(c(a1$x, rev(a1$x)), c(a1$lo, rev(a1$hi)), col = adjustcolor(pal[1], 0.2), border = NA); lines(a1$x, a1$mn, lwd = 2.5, col = pal[1]); points(a1$x, a1$mn, pch = 19, col = pal[1], cex = 0.6)
  abline(v = c(0.85, 0.95), lty = c(3, 2), col = "grey50"); text(0.955, -0.06, expression(rho == 0.95), cex = 0.75, col = "grey40"); text(0.855, -0.06, expression(rho == 0.85), cex = 0.75, col = "grey40")
  emp <- data.frame(rho = c(0.996, 0.993, 0.990, 0.995, 0.960, 0.850, 0.920, 0.730), lb = c("TC24 Mz DTA", "Li18 Sor FT", "Guo20 Rice FT", "Mu22 Sor PH", "Oz25 GY(PHK)", "Wei24 Sor FT", "Han25 Wht FT", "Oz25 GY(PHZ)"), tr = c("FT", "FT", "FT", "PH", "GY", "FT", "FT", "GY"), stringsAsFactors = FALSE)
  epch <- c(FT = 21, PH = 22, GY = 24); ecol <- c(FT = pal[2], PH = pal[3], GY = pal[4]); emp_y <- approx(a1$x, a1$mn, xout = emp$rho, rule = 2)$y
  for (k in seq_len(nrow(emp))) points(emp$rho[k], emp_y[k], pch = epch[emp$tr[k]], bg = ecol[emp$tr[k]], col = "black", cex = 1.6)
  y_off <- c(0.04, -0.04, 0.04, -0.04, 0.04, -0.04, 0.04, -0.04)
  for (k in seq_len(nrow(emp))) text(emp$rho[k], emp_y[k] + y_off[k], emp$lb[k], cex = 0.45, col = "grey30", font = 3)
  legend("bottomleft", legend = c("Simulated (mean, 95% CI)", "Empirical: FT", "PH", "GY"), pch = c(19, 21, 22, 24), pt.bg = c(pal[1], pal[2], pal[3], pal[4]), col = c(pal[1], "black", "black", "black"), cex = 0.7, bty = "n", pt.cex = 1.2); dev.off()
  
  pdf(file.path(outdir, "rho_sweep_power.pdf"), width = 5.5, height = 4.5); par(mar = c(4.8, 5, 1.5, 1), family = "serif", cex.lab = 1.3, cex.axis = 1.1)
  ri <- res_rho[res_rho$n_indep_qtl > 0, ]; if (nrow(ri) > 0) {
    ai <- ag(ri, "rho_target", "power_indep_in_slope"); as2 <- ag(ri, "rho_target", "power_shared_in_slope"); acn <- ag(ri, "rho_target", "power_cn_in_slope"); aap <- ag(ri, "rho_target", "power_ap_in_slope")
    plot(ai$x, ai$mn, type = "n", xlim = c(-0.03, 1.01), ylim = c(0, 1), xlab = expression("Index-mean correlation (" * rho * ")"), ylab = "Detection power (Bonferroni)")
    polygon(c(ai$x, rev(ai$x)), c(ai$lo, rev(ai$hi)), col = adjustcolor(pal[2], 0.15), border = NA); lines(ai$x, ai$mn, lwd = 2.5, col = pal[2]); points(ai$x, ai$mn, pch = 19, col = pal[2], cex = 1.1)
    lines(as2$x, as2$mn, lwd = 2.5, col = pal[1], lty = 2); points(as2$x, as2$mn, pch = 17, col = pal[1], cex = 1.1)
    lines(acn$x, acn$mn, lwd = 2, col = pal[3], lty = 3); points(acn$x, acn$mn, pch = 15, col = pal[3], cex = 1.0)
    lines(aap$x, aap$mn, lwd = 2, col = pal[4], lty = 4); points(aap$x, aap$mn, pch = 18, col = pal[4], cex = 1.2)
    abline(v = c(0.85, 0.95), lty = c(3, 2), col = "grey50")
    legend("topright", legend = c("Independent QTL", "Shared DS QTL", "Conditional neutrality", "Antagonistic pleiotropy"), pch = c(19, 17, 15, 18), col = c(pal[2], pal[1], pal[3], pal[4]), lwd = c(2.5, 2.5, 2, 2), lty = 1:4, cex = 0.65, bty = "n")
  }; dev.off()
  
  pdf(file.path(outdir, "rho_sweep_overlap.pdf"), width = 5.5, height = 4.5); par(mar = c(4.8, 5, 1.5, 1), family = "serif", cex.lab = 1.3, cex.axis = 1.1)
  rs <- res_rho[!is.na(res_rho$overlap_prop), ]; if (nrow(rs) > 0) {
    ao <- ag(rs, "rho_target", "overlap_prop"); plot(ao$x, ao$mn, type = "n", xlim = c(-0.03, 1.01), ylim = c(0, 1), xlab = expression("Index-mean correlation (" * rho * ")"), ylab = "Proportion slope QTL overlapping intercept QTL")
    polygon(c(ao$x, rev(ao$x)), c(ao$lo, rev(ao$hi)), col = adjustcolor(pal[1], 0.2), border = NA); lines(ao$x, ao$mn, lwd = 2.5, col = pal[1]); points(ao$x, ao$mn, pch = 19, col = pal[1], cex = 1.1)
    if (!is.null(res_null)) { ns <- res_null[!is.na(res_null$overlap_prop), ]; if (nrow(ns) > 0) { an <- ag(ns, "rho_target", "overlap_prop"); lines(an$x, an$mn, lwd = 2, col = pal[2], lty = 2); points(an$x, an$mn, pch = 17, col = pal[2], cex = 1.0); legend("bottomright", legend = c("Full architecture", "Null: DS only"), pch = c(19, 17), col = c(pal[1], pal[2]), lwd = 2, lty = 1:2, cex = 0.7, bty = "n") } }
    abline(v = c(0.85, 0.95), lty = c(3, 2), col = "grey50")
  }; dev.off()
}

plot_bootstrap_results <- function(boot_sum, boot_con, res_rho, outdir = ".") {
  key_metrics <- c("power_indep_in_slope", "H2_ratio", "overlap_prop", "cor_logp")
  metric_labs <- c("Power (independent QTL)", expression(H^2 * " ratio (slope/intercept)"), "Slope-intercept QTL overlap", expression("Spearman " * rho * "(-log" * P * ")"))
  if (all(key_metrics %in% boot_sum$metric)) {
    pdf(file.path(outdir, "supp_bootstrap_ci.pdf"), width = 11, height = 8); par(mfrow = c(2, 2), mar = c(4.8, 5, 2.5, 1), family = "serif", cex.lab = 1.1, cex.axis = 1.0)
    for (j in seq_along(key_metrics)) {
      met <- key_metrics[j]; bs <- boot_sum[boot_sum$metric == met, ]; raw <- ag(res_rho, "rho_target", met)
      plot(raw$x, raw$mn, type = "n", xlim = c(-0.03, 1.01), ylim = range(c(raw$lo, raw$hi), na.rm = TRUE), xlab = expression(rho), ylab = metric_labs[j], main = paste0(LETTERS[j], ") ", gsub("_", " ", met)))
      polygon(c(raw$x, rev(raw$x)), c(raw$lo, rev(raw$hi)), col = adjustcolor(pal[6], 0.15), border = NA)
      polygon(c(bs$level, rev(bs$level)), c(bs$ci_lo, rev(bs$ci_hi)), col = adjustcolor(pal[1], 0.3), border = NA)
      lines(bs$level, bs$mean, lwd = 2.5, col = pal[1]); points(bs$level, bs$mean, pch = 19, col = pal[1], cex = 0.5)
      abline(v = c(0.85, 0.95), lty = c(3, 2), col = "grey50")
      if (j == 1) legend("topright", legend = c("Mean (bootstrap)", "95% CI on mean", "95% replicate spread"), col = c(pal[1], adjustcolor(pal[1], 0.3), adjustcolor(pal[6], 0.15)), pch = c(19, 15, 15), pt.cex = c(0.8, 2, 2), cex = 0.6, bty = "n")
    }; dev.off()
  }
}

plot_qg_results <- function(boot_sum, res_rho, outdir = ".") {
  if (all(c("r_predicted", "r_true_genetic", "cor_slope_int") %in% boot_sum$metric)) {
    pdf(file.path(outdir, "supp_theory_overlay.pdf"), width = 6, height = 5); par(mar = c(4.8, 5, 2, 1), family = "serif", cex.lab = 1.3, cex.axis = 1.1)
    bs_emp <- boot_sum[boot_sum$metric == "cor_slope_int", ]; bs_true <- boot_sum[boot_sum$metric == "r_true_genetic", ]; bs_pred <- boot_sum[boot_sum$metric == "r_predicted", ]
    yr <- range(c(bs_emp$ci_lo, bs_emp$ci_hi, bs_true$mean, bs_pred$mean), na.rm = TRUE)
    plot(bs_emp$level, bs_emp$mean, type = "n", xlim = c(-0.03, 1.01), ylim = c(min(yr) - 0.05, max(yr) + 0.05), xlab = expression(rho), ylab = "Correlation (slope, intercept)", main = "Analytical prediction vs. empirical estimation")
    polygon(c(bs_emp$level, rev(bs_emp$level)), c(bs_emp$ci_lo, rev(bs_emp$ci_hi)), col = adjustcolor(pal[1], 0.15), border = NA)
    lines(bs_pred$level, bs_pred$mean, lwd = 2.5, col = pal[3], lty = 1); points(bs_pred$level, bs_pred$mean, pch = 17, col = pal[3], cex = 0.5)
    lines(bs_true$level, bs_true$mean, lwd = 2, col = pal[2], lty = 2); points(bs_true$level, bs_true$mean, pch = 15, col = pal[2], cex = 0.5)
    lines(bs_emp$level, bs_emp$mean, lwd = 2.5, col = pal[1]); points(bs_emp$level, bs_emp$mean, pch = 19, col = pal[1], cex = 0.5)
    abline(v = c(0.85, 0.95), lty = c(3, 2), col = "grey50")
    legend("topleft", legend = c("Empirical", expression("True genetic " * r[G]), "Analytical"), col = c(pal[1], pal[2], pal[3]), pch = c(19, 15, 17), lwd = c(2.5, 2, 2.5), lty = c(1, 2, 1), cex = 0.55, bty = "n"); dev.off()
  }
}

plot_sensitivity_analysis <- function(res_nenv, res_samp, res_marker, res_pop, res_h2 = NULL, res_qtl = NULL, res_hetero = NULL, res_pca = NULL, res_effdist = NULL, res_rho = NULL, res_miss = NULL, res_reps = NULL, res_conv = NULL, res_seed = NULL, res_ldwin = NULL, outdir = ".") {
  sens_panel <- function(res, x_var, y_var, xlab, ylab, main_title, rho_var = "rho_target", log_x = FALSE) {
    rho_levels <- sort(unique(res[[rho_var]])); cols <- if (length(rho_levels) <= 8) pal[seq_along(rho_levels)] else colorRampPalette(c(pal[1], pal[3], pal[2]))(length(rho_levels))
    agg <- do.call(rbind, lapply(split(res, interaction(res[[rho_var]], res[[x_var]])), function(d) { y <- d[[y_var]]; y <- y[!is.na(y)]; data.frame(rho = d[[rho_var]][1], x = d[[x_var]][1], mn = mean(y), se = sd(y) / sqrt(max(length(y), 1)), lo = mean(y) - 1.96 * sd(y) / sqrt(max(length(y), 1)), hi = mean(y) + 1.96 * sd(y) / sqrt(max(length(y), 1)), stringsAsFactors = FALSE) }))
    yr <- range(c(agg$lo, agg$hi), na.rm = TRUE); yr <- yr + c(-1, 1) * diff(yr) * 0.05; if (yr[1] > 0 && yr[1] < 0.1 * diff(yr)) yr[1] <- 0
    plot(NULL, xlim = range(agg$x), ylim = yr, xlab = xlab, ylab = ylab, main = main_title, log = if (log_x) "x" else "")
    for (ri in seq_along(rho_levels)) { sub <- agg[agg$rho == rho_levels[ri], ]; sub <- sub[order(sub$x), ]; if (nrow(sub) < 2) next; polygon(c(sub$x, rev(sub$x)), c(sub$lo, rev(sub$hi)), col = adjustcolor(cols[ri], 0.12), border = NA); lines(sub$x, sub$mn, lwd = 2, col = cols[ri]); points(sub$x, sub$mn, pch = 19, cex = 0.6, col = cols[ri]) }
    legend("topright", legend = paste0("rho = ", rho_levels), col = cols[seq_along(rho_levels)], lwd = 2, pch = 19, cex = 0.5, bty = "n", pt.cex = 0.5)
  }
  
  key_metrics <- c("power_indep_in_slope", "overlap_prop", "H2_ratio", "R_b_gwas_std", "pct_V_b_unique", "cor_slope_int")
  metric_labs <- c("Power (indep)", "Overlap", "H2 ratio", "R_b GWAS", "% unique", "cor(slope, int)")
  
  if (!is.null(res_samp) && nrow(res_samp) > 0) {
    pdf(file.path(outdir, "supp_genotype_sensitivity.pdf"), width = 14, height = 9); par(mfrow = c(2, 3), mar = c(4.5, 5, 2.5, 1), family = "serif", cex.lab = 1.0, cex.axis = 0.9)
    for (mi in seq_along(key_metrics)) if (key_metrics[mi] %in% names(res_samp)) sens_panel(res_samp, "n_geno", key_metrics[mi], "Number of genotypes", metric_labs[mi], paste0(LETTERS[mi], ") ", key_metrics[mi]), log_x = TRUE)
    dev.off()
  }
  if (!is.null(res_pca) && nrow(res_pca) > 0) {
    pdf(file.path(outdir, "supp_pca_lmm_robustness.pdf"), width = 14, height = 9); par(mfrow = c(2, 3), mar = c(4.5, 5, 2.5, 1), family = "serif", cex.lab = 1.0, cex.axis = 0.9)
    for (mi in seq_along(key_metrics)) {
      if (!(key_metrics[mi] %in% names(res_pca))) next
      combos <- unique(paste(res_pca$gwas_method, res_pca$pop_structure)); cols <- pal[seq_along(combos)]
      agg <- do.call(rbind, lapply(split(res_pca, interaction(res_pca$gwas_method, res_pca$pop_structure, res_pca$rho_target)), function(d) { y <- d[[key_metrics[mi]]]; y <- y[!is.na(y)]; data.frame(combo = paste(d$gwas_method[1], d$pop_structure[1]), rho = d$rho_target[1], mn = mean(y), lo = mean(y) - 1.96 * sd(y) / sqrt(max(length(y), 1)), hi = mean(y) + 1.96 * sd(y) / sqrt(max(length(y), 1)), stringsAsFactors = FALSE) }))
      yr <- range(c(agg$lo, agg$hi), na.rm = TRUE); plot(NULL, xlim = c(-0.03, 1.01), ylim = yr, xlab = expression(rho), ylab = metric_labs[mi], main = paste0(LETTERS[mi], ") ", metric_labs[mi]))
      for (ci in seq_along(combos)) { sub <- agg[agg$combo == combos[ci], ]; sub <- sub[order(sub$rho), ]; polygon(c(sub$rho, rev(sub$rho)), c(sub$lo, rev(sub$hi)), col = adjustcolor(cols[ci], 0.12), border = NA); lines(sub$rho, sub$mn, lwd = 2, col = cols[ci], lty = ifelse(grepl("lmm", combos[ci]), 3, ifelse(grepl("pca", combos[ci]), 2, 1))) }
      legend("topright", legend = combos, col = cols, lwd = 2, cex = 0.45, bty = "n", lty = ifelse(grepl("lmm", combos), 3, ifelse(grepl("pca", combos), 2, 1)))
    }
    dev.off()
  }
}

# ============================================================================
# SECTION 14: MAIN EXECUTION
# ============================================================================

main <- function() {
  outdir <- "sim_output_v10"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  env_flag <- function(name, default = TRUE) {
    v <- tryCatch(as.logical(Sys.getenv(name, as.character(default))), error = function(e) default)
    if (length(v) != 1L || is.na(v)) default else v
  }
  RUN_V9  <- env_flag("RN_RUN_V9_SET", TRUE)
  RUN_V10 <- env_flag("RN_RUN_V10_SET", TRUE)
  
  log_file <- file.path(outdir, "run_metadata_and_session.txt")
  cat("========================================================\n", file = log_file)
  cat(" SIMULATION REPRODUCIBILITY LOG (ReactionNormConstraints v10.0.0)\n", file = log_file, append = TRUE)
  cat(" Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)
  cat(sprintf(" SMOKE_TEST = %s; Experiments 1-25 = %s; Experiments 26-29 = %s; workers = %d\n",
              SMOKE, RUN_V9, RUN_V10, N_CORES), file = log_file, append = TRUE)
  cat("========================================================\n\n", file = log_file, append = TRUE)
  cat("--- Default Simulation Parameters ---\n", file = log_file, append = TRUE)
  for (nm in names(default_params)) {
    cat(sprintf("%-15s : %s\n", nm, as.character(default_params[[nm]])), file = log_file, append = TRUE)
  }
  cat("\n--- Session Information ---\n", file = log_file, append = TRUE)
  con <- file(log_file, "a"); writeLines(capture.output(sessionInfo()), con = con); close(con)
  
  cat(strrep("=", 70), "\n")
  cat("  EXECUTION PLAN: Experiments 1-25 (v9.0.0 set) and 26-29 (v10 additions)\n")
  cat(strrep("=", 70), "\n\n")
  t0_global <- proc.time()
  timing <- data.frame(experiment = character(0), seconds = numeric(0), stringsAsFactors = FALSE)
  timed <- function(label, expr) {
    t0 <- proc.time()[3]; val <- expr
    timing[nrow(timing) + 1L, ] <<- list(label, round(proc.time()[3] - t0, 1))
    write.csv(timing, file.path(outdir, "experiment_timing.csv"), row.names = FALSE)
    val
  }

  # Every experiment result is written to CSV so all reported numbers are
  # inspectable from a single run. Filenames are descriptive; see README for the
  # mapping to the experiments and tables of the paper.
  save_csv <- function(obj, name) write.csv(obj, file.path(outdir, name), row.names = FALSE)

  if (RUN_V9) {
  res_rho     <- timed("01 rho_sweep",     run_rho_sweep(rep_mult(1000)));                             save_csv(res_rho,     "rho_sweep.csv")
  res_null    <- timed("02 null_model",    run_null_model(rep_mult(300)));                             save_csv(res_null,    "null_model.csv")
  res_kappa   <- timed("03 kappa",         run_kappa_comparison(rep_mult(200)));                       save_csv(res_kappa,   "kappa_comparison.csv")
  res_ds      <- timed("04 ds_coupling",   run_ds_coupling_sweep(rep_mult(200)));                      save_csv(res_ds,      "ds_coupling_sweep.csv")
  res_envr    <- timed("05 env_range",     run_env_range_sweep(rep_mult(200)));                        save_csv(res_envr,    "env_range_sweep.csv")
  res_nenv    <- timed("06 nenv",          run_nenv_sweep(rep_mult(200)));                             save_csv(res_nenv,    "nenv_sweep.csv")
  res_pop     <- timed("07 pop_structure", run_pop_structure_contrast(rep_mult(150)));                 save_csv(res_pop,     "pop_structure.csv")
  res_arch    <- timed("08 architecture",  run_architecture_comparison(rep_mult(200)));                save_csv(res_arch,    "architecture.csv")
  res_samp    <- timed("09 genotype",      run_sample_size_sweep(rep_mult(100)));                      save_csv(res_samp,    "genotype_sweep.csv")
  res_nl      <- timed("10 nonlinear",     run_nonlinear_test(rep_mult(200)));                         save_csv(res_nl,      "nonlinear_test.csv")
  res_eff     <- timed("11 effect_size",   run_effect_size_sensitivity(rep_mult(200)));                save_csv(res_eff,     "effect_size.csv")
  res_rot     <- timed("12 rotation",      run_rotation_test_validation(rep_mult(50), rep_mult(200))); save_csv(res_rot,     "rotation_test.csv")
  res_multi   <- timed("13 multi_index",   run_multi_index_comparison(rep_mult(150)));                 save_csv(res_multi,   "multi_index.csv")
  res_marker  <- timed("14 marker",        run_marker_density_sweep(rep_mult(100)));                   save_csv(res_marker,  "marker_density.csv")
  res_h2      <- timed("15 heritability",  run_heritability_sweep(rep_mult(100)));                     save_csv(res_h2,      "heritability_sweep.csv")
  res_qtl     <- timed("16 qtl_count",     run_qtl_count_sweep(rep_mult(100)));                        save_csv(res_qtl,     "qtl_count.csv")
  res_hetero  <- timed("17 hetero_res",    run_heterogeneous_residuals(rep_mult(100)));                save_csv(res_hetero,  "heterogeneous_residuals.csv")
  res_pca     <- timed("18 pca_lmm",       run_pca_gwas_robustness(rep_mult(100)));                    save_csv(res_pca,     "pca_lmm_gwas.csv")
  res_effdist <- timed("19 effect_dist",   run_effect_distribution(rep_mult(100)));                    save_csv(res_effdist, "effect_distribution.csv")
  res_miss    <- timed("20 missing",       run_missing_data_sweep(rep_mult(100)));                     save_csv(res_miss,    "missing_data.csv")
  res_reps    <- timed("21 replication",   run_replication_sweep(rep_mult(100)));                      save_csv(res_reps,    "replication_sweep.csv")
  res_conv    <- timed("22 convergence",   run_convergence_diagnostic(if (SMOKE) 0.02 else 1));        save_csv(res_conv,    "convergence.csv")
  res_seed    <- timed("23 seed",          run_seed_sensitivity(rep_mult(100)));                       save_csv(res_seed,    "seed_sensitivity.csv")
  res_ldwin   <- timed("24 ld_window",     run_ld_window_sweep(rep_mult(100)));                        save_csv(res_ldwin,   "ld_window.csv")
  res_maf     <- timed("25 maf",           run_maf_sweep(rep_mult(100)));                              save_csv(res_maf,     "maf_sweep.csv")
  }
  
  if (RUN_V10) {
  res_so    <- timed("26 slope_only_qtl", run_slope_only_qtl_sweep(rep_mult(200)));  save_csv(res_so,   "slope_only_qtl.csv")
  save_csv(summarize_v10(res_so, c("kappa", "so_axis", "so_center", "rho_target")), "summary_slope_only_qtl.csv")
  res_core  <- timed("27 core_simulation", run_core_simulation(rep_mult(200)));     save_csv(res_core, "core_simulation.csv")
  save_csv(summarize_v10(res_core, "rho_target"), "summary_core_simulation.csv")
  res_wev   <- timed("28 within_env_variance", run_within_env_variance(rep_mult(100))); save_csv(res_wev, "within_env_variance.csv")
  save_csv(summarize_v10(res_wev, c("part", "gxe_mean_axis", "sigma_res", "n_reps", "var_a_mean",
                                    "env_range", "n_env", "rho_target")), "summary_within_env_variance.csv")
  res_pve   <- timed("29 locus_pve", run_locus_pve(rep_mult(200)));                 save_csv(res_pve,  "locus_pve.csv")
  save_csv(summarize_v10(res_pve, c("n_so_qtl", "rho_target")), "summary_locus_pve.csv")
  }

  if (RUN_V9) {
  cat("\n  Generating summaries...\n")
  rho_sum <- summarize_rho_sweep(res_rho);                                  save_csv(rho_sum,     "summary_rho_sweep.csv")
  pca_sum <- summarize_by(res_pca, c("rho_target", "gwas_method", "pop_structure")); save_csv(pca_sum, "summary_pca_lmm_gwas.csv")

  n_boot <- rep_mult(2000)
  cat(sprintf("\n  Bootstrap inference (n_boot = %d)...\n\n", n_boot))
  set.seed(20260224L)
  boot_rho <- bootstrap_summary(res_rho, "rho_target", n_boot = n_boot)
  boot_con <- bootstrap_contrasts(res_rho, "rho_target", n_boot = n_boot)
  
  cat("\n  Generating publication figures...\n")
  plot_results(res_rho, res_null, res_arch, res_ds, res_nenv, res_pop, res_rot, res_multi, outdir)
  plot_bootstrap_results(boot_rho, boot_con, res_rho, outdir)
  plot_qg_results(boot_rho, res_rho, outdir)
  plot_sensitivity_analysis(res_nenv, res_samp, res_marker, res_pop, res_h2, res_qtl, res_hetero, res_pca, res_effdist, res_rho, res_miss, res_reps, res_conv, res_seed, res_ldwin, outdir)
  }
  
  elapsed <- (proc.time() - t0_global)[3]
  cat(sprintf("\n%s\n  COMPLETE. Total time: %.1f min (%.0f sec)\n", strrep("=", 70), elapsed / 60, elapsed))
  cat(sprintf("  Output directory: %s\n", normalizePath(outdir)))
}

# ============================================================================
# SECTION 15: STANDALONE ROTATION TEST FOR SUPPLEMENT
# ============================================================================

standalone_rotation_test <- function(Y, X, x, n_perm = 500, alpha_bonf = 0.05) {
  stopifnot(is.matrix(Y), is.matrix(X), length(x) == ncol(Y))
  n  <- nrow(Y); nm <- ncol(X); th <- alpha_bonf / nm
  
  xc <- x - mean(x); ym <- rowMeans(Y, na.rm = TRUE); vx <- sum(xc^2)
  slopes <- as.vector((Y - ym) %*% xc) / vx; intercepts <- ym - slopes * mean(x)
  rho_obs <- cor(x, colMeans(Y, na.rm = TRUE))
  
  yc_int <- intercepts - mean(intercepts); yc_slp <- slopes - mean(slopes)
  ssx <- colSums(X^2); mono <- ssx < 1e-10; ssx[mono] <- 1
  
  calc_pv <- function(yc) {
    beta <- as.vector(crossprod(X, yc)) / ssx
    ssr <- pmax(sum(yc^2) - beta^2 * ssx, 1e-10)
    se <- sqrt((ssr / (n - 2)) / ssx)
    pv <- 2 * pt(abs(beta / se), df = n - 2, lower.tail = FALSE)
    pv[mono] <- 1; pv
  }
  
  si_obs <- which(calc_pv(yc_int) < th); ss_obs <- which(calc_pv(yc_slp) < th)
  unique_obs <- length(setdiff(ss_obs, si_obs))
  
  noise_sd <- sqrt(max(1 - rho_obs^2, 1e-10))
  Y_null <- matrix(rnorm(n * n_perm, 0, noise_sd), n, n_perm)
  Y_null <- sweep(Y_null, 1, scale(intercepts, scale = FALSE)[, 1] * rho_obs, "+")
  Y_null_c <- scale(Y_null, center = TRUE, scale = FALSE)
  
  ssy <- colSums(Y_null_c^2); XY <- crossprod(X, Y_null_c)
  beta <- XY / ssx; ssr <- sweep(-(beta * XY), 2, ssy, "+"); ssr[ssr < 1e-10] <- 1e-10
  se <- sqrt(sweep(ssr / (n - 2), 1, ssx, "/"))
  pv_mat <- 2 * pt(abs(beta / se), df = n - 2, lower.tail = FALSE); pv_mat[mono, ] <- 1
  
  if (length(si_obs) > 0) null_unique <- colSums(pv_mat[-si_obs, , drop = FALSE] < th) else null_unique <- colSums(pv_mat < th)
  
  list(rho_observed = rho_obs, n_unique_slope = unique_obs, null_mean = mean(null_unique), null_q95 = unname(quantile(null_unique, 0.95)), pvalue = (sum(null_unique >= unique_obs) + 1) / (n_perm + 1))
}

if (!interactive()) main()

