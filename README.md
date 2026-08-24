# Causal Fairness for Survival Analysis

This repository contains the code for the paper *Causal Fairness for Survival Analysis*.

1. Code for reproducing the experiments is in `scripts/anzics-survival.R`.
2. Model-based estimation is implemented in `fair_surv()` in `r/fair-surv.R`.
3. Doubly robust estimation is implemented in `one_step_debias_surv()` in `r/one-step-debias-surv.R`.
4. Synthetic data generators (log-normal and Weibull, dispatched via S3 on `class(dgm)`) are in `r/synth-shared.r`, `r/synth-lognormal.r`, and `r/synth-weibull.R`.
5. Synthetic validation experiments (coverage, competing risks, misspecification robustness) are in `scripts/nic-synth.R`, `scripts/cr-synth.R`, and `scripts/nic-misspec.R`.