# Landmark-Based Dynamic Prediction of Mean Residual Life (MRL)

Dynamic prediction of **mean residual lifetime (MRL)** using a landmark supermodel framework, comparing **pseudo-observation (PO)** and **inverse probability of censoring weighting (IPCW)** approaches. The method is illustrated on the `pbc2` primary biliary cirrhosis dataset.

## Overview

Mean residual lifetime at a prediction time *s* is the expected remaining survival time given that a subject is still alive at *s*, restricted to a horizon τ:

> MRL(s) = E[ min(T, τ) − s | T > s ]

This repository implements and compares two estimation strategies for MRL within a landmark framework, and contrasts a single **dynamic supermodel** (one model fitted across stacked landmark datasets, with time-varying effects via landmark interactions) against a series of **static models** (one model per landmark time).

- **Pseudo-observation (PO) approach** — jackknife pseudo-values of the restricted residual lifetime, regressed via GEE.
- **IPCW approach** — observed residual times reweighted by the inverse of the estimated censoring survival probability.

Both approaches are wrapped in a landmark supermodel that allows covariate effects to vary smoothly with prediction time through linear and quadratic landmark terms.

## Method summary

1. **Data preparation.** Apply a finite horizon τ, and define the censoring/event indicators used by the IPCW weights.
2. **Landmark super dataset.** For a grid of landmark times, subjects still at risk are stacked into a long-format super dataset (via `cutLM`), with the residual time and pseudo-observations computed within each landmark stratum.
3. **Time-varying effects.** Covariates are interacted with scaled landmark polynomials (`.t0`, `.t1`, `.t2`) so effects can evolve over prediction time.
4. **Supermodels.** A dynamic PO model and a dynamic IPCW model are fitted with GEE (`geeglm`, independence working correlation). Static per-landmark models are fitted for comparison.
5. **Prediction & visualization.** Coefficient trajectories over landmark time, and patient-level MRL predictions (static vs. dynamic) with confidence bands.
6. **Validation.** Monte-Carlo train/test splits evaluate discrimination (a residual-life C-index) and an IPCW-weighted absolute prediction error across landmark times.

## Requirements

- R packages:
  - [`JM`](https://cran.r-project.org/package=JM) — provides the `pbc2` dataset
  - [`pseudo`](https://cran.r-project.org/package=pseudo) — pseudo-observation calculation (`pseudomean`)
  - [`dynpred`](https://cran.r-project.org/package=dynpred) — landmarking utilities (`cutLM`)
  - [`gee`](https://cran.r-project.org/package=gee) and [`geepack`](https://cran.r-project.org/package=geepack) — GEE models (`gee`, `geeglm`)
  - [`survival`](https://cran.r-project.org/package=survival) — `survfit` for the IPCW weights

Key settings you may want to adjust near the top of the script:

| Setting | Meaning |
|---|---|
| `tau` | Restriction horizon for residual lifetime |
| `sL` | Maximum landmark (prediction) time considered |
| `LMs` | Grid of landmark times |
| `fixed` | Baseline (time-fixed) covariates |
| `varying` | Time-varying covariates |
| `u` | Number of Monte-Carlo validation replicates |

## Outputs

- Coefficient summaries (estimate, SE, Z, p-value) for the dynamic and static PO/IPCW models.
- Coefficient-evolution plots showing how covariate effects change over prediction time.
- Patient-level MRL prediction plots comparing the static fixed models with the dynamic supermodel, including confidence intervals.
- Validation plots of the C-index and prediction error versus prediction time for all four model variants.
