# ReactionNormConstraints

Companion code for:

**Detection power for plasticity loci when the environmental index is chosen to match environment means**
Shawn D. Jenkins and George L. Graef.

Version 10.0.0. See [CHANGES.md](CHANGES.md).

## What this is

CERIS-JGRA regresses each genotype's performance on a climate index chosen, among many candidate climate variables and growth windows, for its correlation (rho) with the environment means, and then maps the fitted slopes and intercepts by GWAS. Any such index can be written as a part along the environment means, with weight rho, plus a part along other environmental variation, with weight tau, where tau cannot exceed sqrt(1 - rho^2). This repository holds the simulations, checks and reanalysis behind the paper: how slope-GWAS power for different classes of loci changes with rho, an analytic prediction of that power, and a reimplementation of the CERIS search on the published sorghum data.

The simulation is one script, `sim_reaction_norm_gwas_v10.R`. Version 10.0.0 adds a slope-only QTL class with zero intercept (or zero mean) effect, an optional simplified simulation (independent loci, fixed allele frequency, fixed effect sizes), an analytic power prediction that includes within-environment genotypic and residual variance, per-locus variance explained, and Experiments 26 to 29. All additions are off by default, and Experiments 1 to 25 reproduce v9.0.0 replicate for replicate.

| File | Purpose |
|---|---|
| `sim_reaction_norm_gwas_v10.R` | Simulation framework, Experiments 1 to 29 |
| `reaction_norm_audit.R` | Additional runs in four stages: diagnostics, the five-architecture comparison, a finer coupling grid, and simulated against published -log10(P) correlations |
| `run_full_v10.R`, `run_full_v10.ps1` | Driver and PowerShell launcher for the full reruns |
| `tools/sorghum_ceris_scan.R` | CERIS search on the Li et al. (2018) sorghum data: 27,144 candidates, rho and largest possible tau |
| `tools/closed_form_design_experiments.R` | Closed-form power at the 126 design points of the kappa, panel-size and environment-number experiments, beside observed power |
| `tools/derivation_check.R` | Numerical check of the power derivation |
| `tools/reproduce_v9_check.R` | Reproduction check against v9.0.0 output |
| `tools/benchmark_runtime.R` | Per-replicate cost and runtime estimate for your machine |
| `sim_output_v10/`, `audit_output_v10/`, `sorghum_scan/` | Output of the v10 reference run |
| `sim_output_v9/` | v9.0.0 summary data (unchanged) |

## Requirements

- R >= 3.6.0
- Base R only. The `parallel` package ships with R, so there is nothing to install.

## Reproduce

Run from the repository root.

Everything (simulation, additional runs, derivation check), on Windows PowerShell:

```powershell
.\run_full_v10.ps1              # full scale
.\run_full_v10.ps1 -Smoke       # pipeline check, replicate counts scaled to ~2%
.\run_full_v10.ps1 -Parts "sim_v10,audit" -Cores 24
.\run_full_v10.ps1 -DryRun      # plan and runtime estimate only
```

or with any shell:

```sh
Rscript run_full_v10.R [--parts=sim_v9,sim_v10,audit,check] [--smoke] [--cores=K] [--dry-run]
```

The simulation alone:

```sh
Rscript sim_reaction_norm_gwas_v10.R
SMOKE_TEST=TRUE Rscript sim_reaction_norm_gwas_v10.R
RN_RUN_V9_SET=FALSE Rscript sim_reaction_norm_gwas_v10.R     # Experiments 26-29 only
```

The reference run took 138 minutes with 30 workers on a 32-thread workstation (R 4.6.1): 124 minutes for Experiments 1 to 29 (Experiment 1 alone 41 minutes), 11 minutes for the additional runs and 4 minutes for the derivation check. `tools/benchmark_runtime.R` estimates these for another machine. A smoke run of everything takes roughly 15 to 20 minutes. Simulation output is written to `./sim_output_v10/`, output of the additional runs to `./audit_output_v10/`, logs to `./logs/`.

## Reproducibility

- The RNG regime is pinned (`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`). This is R's default from 3.6.0 onward, so pinning does not change the published numbers; it makes them stable across R installations at that version or later.
- Per-replicate seeds are deterministic, so results do not depend on the number of cores.
- `sim_output_v10/run_metadata_and_session.txt` records the default parameters and full `sessionInfo()` for the run, and `sim_output_v10/experiment_timing.csv` the wall time of each experiment.
- Across R builds, detection outcomes (power, counts, overlap) are identical and continuous summaries agree to about 1e-9 (reference BLAS against OpenBLAS).

## Output

A single run produces, in `./sim_output_v10/`:

- 29 experiment result tables (CSV)
- 6 summary tables (`summary_rho_sweep.csv`, `summary_pca_lmm_gwas.csv`, `summary_slope_only_qtl.csv`, `summary_core_simulation.csv`, `summary_within_env_variance.csv`, `summary_locus_pve.csv`)
- 7 diagnostic figure PDFs
- `run_metadata_and_session.txt`

`sim_output_v9/` holds the v9.0.0 summary data and is kept unchanged. The v10 reference run is committed in `sim_output_v10/` and `audit_output_v10/`. The largest table, `rho_sweep.csv` (81 MB), is committed compressed as `rho_sweep.csv.gz`; read it with `read.csv(gzfile("sim_output_v10/rho_sweep.csv.gz"))`.

### Experiments

Numbered in run order, as in Table S1 of the paper.

| # | CSV | What it varies |
|---|-----|----------------|
| 1 | `rho_sweep.csv` | rho over a fine grid (0 to 0.999), 1000 reps (2000 at 0.94 and 0.95). Primary result. |
| 2 | `null_model.csv` | No QTL of any class, rho fine grid. Calibration and empirical false-positive rate. |
| 3 | `kappa_comparison.csv` | rho (coarse) x kappa in {0, 0.25, 0.5, 0.75, 1} (7 x 5 = 35 cells). |
| 4 | `ds_coupling_sweep.csv` | rho (coarse) x differential-sensitivity coupling in {0, 0.3, 0.6, 0.9, 1}. |
| 5 | `env_range_sweep.csv` | rho (coarse) x environmental range {wide, narrow}. |
| 6 | `nenv_sweep.csv` | rho (coarse) x number of environments in {3, 4, 5, 6, 8, 12, 16, 20}. |
| 7 | `pop_structure.csv` | rho (coarse) x population {diverse, biparental, nam}. |
| 8 | `architecture.csv` | rho (coarse) across QTL-architecture configurations. |
| 9 | `genotype_sweep.csv` | rho (sensitivity) x n_geno in {50, 100, 200, 500, 1000, 2000, 5000}. |
| 10 | `nonlinear_test.csv` | rho (coarse) x nonlinear reaction norm {FALSE, TRUE}. |
| 11 | `effect_size.csv` | rho in {0, 0.5, 0.85, 0.95, 0.996} x independent-QTL effect variance in {0.01, 0.02, 0.04, 0.08}. |
| 12 | `rotation_test.csv` | Rotation-based permutation test over rho (coarse): slope-unique discovery against the null. |
| 13 | `multi_index.csv` | Multi-index GWAS comparison over rho (coarse). |
| 14 | `marker_density.csv` | rho (sensitivity) x n_markers in {500, 1000, 2000, 5000, 10000}. |
| 15 | `heritability_sweep.csv` | rho (sensitivity) x residual SD in {0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 8.0}. |
| 16 | `qtl_count.csv` | rho (sensitivity) across total-QTL configurations (polygenicity, ~10 to ~1000 QTL). |
| 17 | `heterogeneous_residuals.csv` | rho (sensitivity) x homoscedastic vs heteroscedastic residuals. |
| 18 | `pca_lmm_gwas.csv` | rho (sensitivity) x method {simple, pca, lmm} x population {diverse, nam}. |
| 19 | `effect_distribution.csv` | Equal vs unequal effect-size distribution over rho (sensitivity). |
| 20 | `missing_data.csv` | rho (sensitivity) x missing-data fraction in {0, 0.05, 0.1, 0.2}. |
| 21 | `replication_sweep.csv` | rho (sensitivity) x within-environment reps in {1, 2, 3, 4}. |
| 22 | `convergence.csv` | Monte Carlo convergence at rho in {0.5, 0.85, 0.996} across increasing replicate counts. |
| 23 | `seed_sensitivity.csv` | Five seed groups x rho in {0.5, 0.85, 0.996}. |
| 24 | `ld_window.csv` | rho (sensitivity) x LD detection window in {0, 2, 5, 10}. |
| 25 | `maf_sweep.csv` | rho (sensitivity) x minor-allele-frequency threshold in {0.01, 0.05, 0.10, 0.20, 0.35}. |
| 26 | `slope_only_qtl.csv` | Default architecture plus 5 slope-only loci: rho (coarse) x axis {index, driver}; centering {intercept, mean} at kappa = 0.5. |
| 27 | `core_simulation.csv` | Simplified core simulation (independent loci, allele frequency 0.30, fixed effects) plus 5 index-axis slope-only loci, rho (coarse). |
| 28 | `within_env_variance.csv` | Core simulation: residual SD x replicates x mean-axis GxE; main-effect variance x environmental range x number of environments; rho (sensitivity). Observed and predicted power. |
| 29 | `locus_pve.csv` | Per-locus variance explained by class (per environment, genotype means, slope, intercept), with and without slope-only loci. |

### Figures

The script writes seven diagnostic PDFs summarizing the primary rho sweep and the robustness checks, named by content. The figures in the paper (Figures 1 to 4 and S1 to S12) were drawn separately from the output tables.

| PDF | Content |
|-----|---------|
| `rho_sweep_cor.pdf` | Slope-intercept effect-size correlation vs rho. |
| `rho_sweep_power.pdf` | Power for independent-axis and differential-sensitivity loci vs rho, with published studies overlaid. |
| `rho_sweep_overlap.pdf` | Slope-intercept QTL overlap proportion vs rho. |
| `supp_bootstrap_ci.pdf` | Bootstrap confidence intervals (2x2). |
| `supp_theory_overlay.pdf` | Analytical prediction vs empirical estimate of the slope-intercept correlation. |
| `supp_genotype_sensitivity.pdf` | Genotype and sample-size sensitivity panels (2x3). |
| `supp_pca_lmm_robustness.pdf` | PCA and LMM robustness (2x3). |

## Sorghum reanalysis

The sorghum data are public: https://github.com/jmyu/CERIS-JGRA-Workshop (folder `1Sorghum`, commit 730af93). Clone that repository, then from this repository's root:

```sh
Rscript tools/sorghum_ceris_scan.R path/to/CERIS-JGRA-Workshop/1Sorghum sorghum_scan
```

The trait is flowering time in growing degree days (`FTgdd`), averaged over the 237 lines in each of seven environments. The search covers day length, growing degree days, photothermal time and photothermal ratio over every window of at least seven days within 122 days after planting (4 x 6,786 = 27,144 candidates). It runs in about a second and writes `sorghum_scan/sorghum_ceris_scan.csv` and a summary.

## Meta-analysis

`cross_study.csv` lists the 47 study-trait combinations compiled from 27 publications. Its values are taken from the cited papers; Table S2 of the paper gives the same rows with verified DOIs.

## Citation

Jenkins, S.D., and G.L. Graef. Detection power for plasticity loci when the environmental index is chosen to match environment means. Submitted to G3: Genes|Genomes|Genetics.

Code: ReactionNormConstraints v10.0.0, DOI to be added when the release is archived on Zenodo. v9.0.0 is archived at Zenodo, DOI 10.5281/zenodo.21610339.

## License

MIT. See [LICENSE](LICENSE).
