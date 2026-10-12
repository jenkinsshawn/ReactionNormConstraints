# Changes

## 10.0.0

This release adds analyses and options. With every new option at its default,
Experiments 1-25 reproduce v9.0.0 replicate for replicate: the same seeds give the
same values in every v9.0.0 column. The full-scale Experiment 1 summary matches
v9.0.0 to 1e-12.

### Added

- **No-intercept slope QTL class** (`n_so_qtl`, `var_a_so`, `so_axis`, `so_center`;
  Section 4B). A slope-only locus adds `g_i * v_j` to genotype i in environment j.
  `so_axis = "index"` builds `v` from the index `x`; `so_axis = "driver"` builds it
  from the environmental driver `h`. `so_center = "intercept"` gives exactly zero
  effect on the intercept at `x = 0`; `so_center = "mean"` gives exactly zero effect
  on the environment-averaged genotype mean. The two coincide when `mean(x) = 0`,
  which holds at `kappa = 1` up to floating-point rounding. Positions and effects
  come from an independent RNG substream, so switching the class on leaves every
  other locus, effect, environment, residual and index draw of a replicate
  unchanged. Off by default (`n_so_qtl = 0`).
- **Simplified core simulation** (`core_sim`, `core_maf`, `fixed_effects`;
  Section 2B): independent loci, one allele frequency for every marker (default
  0.30), fixed effect magnitudes with random sign, no LD. Off by default.
- **Analytic power prediction** (`predict_power`; Section 7C):
  `analytic_power_conditional()` gives per-locus non-centrality and power from the
  slope and intercept contrasts of the realized index and the true generating
  values; `closed_form_power()` gives the expectation over environments. See
  Appendix A of the paper.
- **Per-locus variance explained by QTL class** (`pve_cols`; `locus_pve()`), per
  environment, on genotype means, on the estimated slope and on the estimated
  intercept.
- **Experiments 26-29**: slope-only QTL class (26), core simulation (27),
  within-environment variance (28), per-locus variance explained (29), with
  summary tables.
- **Sorghum CERIS search** (`tools/sorghum_ceris_scan.R`): the CERIS search on the
  public Li et al. (2018) sorghum data (github.com/jmyu/CERIS-JGRA-Workshop, folder
  `1Sorghum`); writes rho and the largest possible tau for all 27,144 candidates to
  `sorghum_scan/`.
- **Closed-form power at the design points** (`tools/closed_form_design_experiments.R`):
  closed-form power for independent-axis loci at the 126 cells of the kappa,
  panel-size and environment-number experiments, beside observed power.
- `reaction_norm_audit.R`: additional runs in four stages (diagnostics, five QTL
  architectures, a finer coupling grid, simulated against published -log10(P)
  correlations), written to `./audit_output_v10/`. Its extra columns come from
  `run_replicate()` when `audit_cols = TRUE`, including `H2_slope_x`, slope
  heritability scored against the exact genotypic slope on x.
- `tools/derivation_check.R`: numerical check of the power derivation.
- `tools/reproduce_v9_check.R`: re-runs a committed v9.0.0 summary row at full
  scale and seed-matched v9.0.0 replicates, and compares every column.
- `tools/benchmark_runtime.R`: per-replicate cost under load and a full-run
  runtime estimate for the local machine.
- `run_full_v10.R` and `run_full_v10.ps1`: driver and PowerShell launcher for the
  full reruns, with per-part logs and runtime estimates.
- Run control by environment variable: `RN_RUN_V9_SET`, `RN_RUN_V10_SET`,
  `RN_N_CORES`. `sim_output_v10/experiment_timing.csv` records per-experiment
  wall time.
- Reference run: full-scale output (R 4.6.1, 30 workers) in `sim_output_v10/`,
  `audit_output_v10/` and `sorghum_scan/`; `rho_sweep.csv` is committed
  gzip-compressed.

### Changed

- `sim_reaction_norm_gwas_v9.R` is renamed `sim_reaction_norm_gwas_v10.R` and
  writes to `./sim_output_v10/`. `sim_output_v9/` is kept unchanged as the v9.0.0 data.
- `run_experiment()` attaches the design factors that `run_replicate()` does not
  record (for example `kappa`, `env_range`, `gwas_method`, `sigma_res`, `n_reps`,
  `n_markers`, `maf_min`, `pct_missing`, `hetero_res`, `nonlinear`,
  `var_a_indep`) as extra columns, so each experiment table can be grouped by the
  factor it varies. `summary_pca_lmm_gwas.csv` now reports the simple, PCA and
  mixed-model tests separately; the v9.0.0 file pools them. Values in existing
  columns are unchanged.
- Optional v10 columns that are NA in every row are dropped, so default
  experiment tables keep the v9.0.0 column set apart from the design columns
  above.
- `run_experiment()` handles experiments with a single task, which `SMOKE_TEST`
  mode produces for Experiments 22-23. Chunking does not affect results.
- The worker task function no longer carries the `run_experiment()` frame when
  it is shipped to workers.
- `generate_phenotypes()` also returns the per-environment residual SD `sig_j`
  (no change to draws).
- `.gitignore` and `.gitattributes` replace the undotted `gitignore` and
  `gitattributes`, so Git applies them.

### Not changed

- Phenotype model, index construction, reaction norm estimation, GWAS, metrics,
  seeds and RNG regime for Experiments 1-25.
