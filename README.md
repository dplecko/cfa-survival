# Causal Fairness for Survival Analysis

This repository contains the code for the paper *Causal Fairness for Survival Analysis*.

## Methods (`r/`)

- Doubly robust (one-step) estimation: `one_step_debias_surv()` in `r/one-step-debias-surv.R`, including the copula-graphic sensitivity analysis for informative censoring (envelope, Routes I and II, with partial-identification CIs).
- Model-based estimation: `fair_surv()` in `r/fair-surv.R`.
- Synthetic data generators (log-normal and Weibull, dispatched via S3 on `class(dgm)`): `r/synth-shared.r`, `r/synth-lognormal.r`, `r/synth-weibull.R`.

## Experiments (`scripts/`)

| Paper result | Script |
|---|---|
| ANZICS death analysis (main text) | `scripts/anzics-survival.R`, task 1 |
| ANZICS readmission with informative-censoring sensitivity | `scripts/anzics-survival.R`, task 2 |
| ANZICS competing risks (readmission and death) | `scripts/anzics-survival.R`, task 2 |
| Seed robustness of the death analysis | `scripts/anzics-survival.R`, tasks 1 and 3-11, then `scripts/anzics-seeds.R` |
| Z/W partition sensitivity | `scripts/anzics-survival.R`, tasks 1, 12, 13, then `scripts/anzics-zw-sensitivity.R` |
| Coverage, non-informative censoring | `scripts/nic-synth.R` |
| Coverage, competing risks | `scripts/cr-synth.R` |
| Multiple robustness (RMSE vs n) | `scripts/nic-misspec.R` |
| Further misspecification arms (S2, C1a, C1b, C2) | `scripts/misspec-synth.R` (tasks, then `agg`) |
| Informative censoring with known truth (Routes I and II) | `scripts/ic-synth.R` (tasks, then `agg`) |

`scripts/anzics-survival.R`, `scripts/ic-synth.R` and `scripts/misspec-synth.R` are task arrays: each task is selected by `SLURM_ARRAY_TASK_ID` (or `SIM_TASK` for the synthetic scripts); `SIM_LIST=1` prints the task table. Without a task id, `scripts/anzics-survival.R` runs task 1 on a 5,000-patient subsample. The ANZICS data are not included.
