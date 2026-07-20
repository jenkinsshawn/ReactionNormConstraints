# ReactionNormConstraints

Companion code for:

**When Does Reaction Norm GWAS Discover Plasticity? The Residual Channel in Environmental Index-Based Genetic Dissection**
Shawn D. Jenkins and George L. Graef. Submitted to *GENETICS* (GSA).

## What this is

A simulation framework that tests when a reaction-norm GWAS built on an environmental index can separate phenotypic plasticity (reaction-norm slope) from the trait mean (intercept), as a function of the index-mean correlation rho. The central result is algebraic: maximizing rho, which is the objective of index-construction methods such as CERIS, caps the residual-channel loading at tau = sqrt(1 - rho^2), so a plasticity-specific locus retains a vanishing fraction of its non-centrality as rho approaches 1.

Everything is one script: `sim_reaction_norm_gwas_v9.R`.

## Requirements

- R >= 3.6.0
- Base R only. The `parallel` package ships with R, so there is nothing to install.

## Reproduce

Run from the repository root.

Full sweep:

```sh
Rscript sim_reaction_norm_gwas_v9.R
```

Quick pipeline check (scales every replicate, permutation, and bootstrap count to ~2%, finishes in minutes):

```sh
SMOKE_TEST=TRUE Rscript sim_reaction_norm_gwas_v9.R
```

On Windows PowerShell:

```powershell
$env:SMOKE_TEST = "TRUE"; Rscript sim_reaction_norm_gwas_v9.R
```

Full-sweep runtime is roughly 2 to 3 hours on a 12-core / 128 GB machine. Output is written to a `./sim_output_v9/` directory created at run time.

## Reproducibility

- The RNG regime is pinned (`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`). This is R's default from 3.6.0 onward, so pinning does not change the published numbers; it makes them stable across R installations at that version or later.
- Per-replicate seeds are deterministic, so results do not depend on the number of cores.
- `sim_output_v9/run_metadata_and_session.txt` records the default parameters and full `sessionInfo()` for the run.

## Output

A single run produces, in `./sim_output_v9/`:

- 25 experiment result tables (CSV)
- 2 summary tables (`summary_rho_sweep.csv`, `summary_pca_lmm_gwas.csv`)
- 7 diagnostic figure PDFs
- `run_metadata_and_session.txt`

The CSVs in this directory are the summary data files the manuscript's Data Availability statement points to. Run the script once and commit `sim_output_v9/` so the data ships with the code. With the RNG pinned, re-running reproduces byte-identical files, so the committed copy stays clean.

### Experiments

Numbered in run order, which follows the manuscript's experiment table.

| # | CSV | What it varies |
|---|-----|----------------|
| 1 | `rho_sweep.csv` | rho over a fine grid (0 to 0.999), 1000 reps. Primary result. |
| 2 | `null_model.csv` | No QTL of any class, rho fine grid. Calibration and empirical false-positive rate. |
| 3 | `kappa_comparison.csv` | rho (coarse) x kappa in {0, 0.25, 0.5, 0.75, 1} (7 x 5 = 35 cells). |
| 4 | `ds_coupling_sweep.csv` | rho (coarse) x directional-selection coupling in {0, 0.3, 0.6, 0.9, 1}. |
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

### Figures

The script writes seven diagnostic PDFs summarizing the primary rho sweep and the robustness checks. These are the simulation's own renders, named by content. The manuscript's numbered figures (main text Figures 1 to 5 and the supplemental figures) are produced by separate figure-generation scripts and are not the files below.

| PDF | Content |
|-----|---------|
| `rho_sweep_cor.pdf` | Slope-intercept effect-size correlation vs rho. |
| `rho_sweep_power.pdf` | Plasticity-specific and differentially-sensitive power vs rho, with published studies overlaid. |
| `rho_sweep_overlap.pdf` | Slope-intercept QTL overlap proportion vs rho. |
| `supp_bootstrap_ci.pdf` | Bootstrap confidence intervals (2x2). |
| `supp_theory_overlay.pdf` | Analytical prediction vs empirical estimate of the slope-intercept correlation. |
| `supp_genotype_sensitivity.pdf` | Genotype and sample-size sensitivity panels (2x3). |
| `supp_pca_lmm_robustness.pdf` | PCA and LMM robustness (2x3). |

## Scope and empirical analyses

This repository contains the simulation framework. The sorghum reanalysis and the cross-study meta-analysis reported in the paper draw on data from the cited studies, not on anything generated here. The empirical points overlaid in Figure 2b are taken from published work (Tibbs-Cortes et al. 2024; Li et al. 2018; Guo et al. 2020; and the additional studies listed in the manuscript) and are not simulation output.

## Citation

Jenkins, S.D., and G.L. Graef. 2026. When Does Reaction Norm GWAS Discover Plasticity? The Residual Channel in Environmental Index-Based Genetic Dissection. *GENETICS*. DOI: to be minted upon acceptance.

## License

MIT. See [LICENSE](LICENSE).
