
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

DGM_SEED <- 2026

# settings
pz <- 2
q <- 2
k <- 2
n <- 10000
nseed <- 96
n_truth <- 1e6
effects <- c("tv", "ctfde", "ctfie", "ctfse")
events <- 1:2

X <- "majority"
Z <- paste0("z", seq_len(pz))
W <- paste0("w", seq_len(q))
event_var <- "event"
time_var <- "event_time"

# fix one DGP
set.seed(DGM_SEED)
T_par  <- list(mu0 = -0.15, muX = -0.15,
               muZ = c(-0.08, 0.07),
               muW = c(-0.11, -0.15), shape = 1.5, scale = 1)

T2_par <- list(mu0 = -0.15, muX =  0.15,
               muZ = c( 0.07, -0.06),
               muW = c( 0.13,  0.13), shape = 1.5, scale = 1)

C_par  <- list(mu0 = -0.85, muX = -0.08,
               muZ = c(-0.04, 0.04),
               muW = c(-0.03, 0.03), shape = 1.2, scale = 1)
g0 <- gen_surv_weibull(n = 1, k = k, pz = pz, q = q, C_par = C_par,
               T_par = T_par, T2_par = T2_par, seed = DGM_SEED,
               pW = 0.75, alpha = rep(0.75, q))
par <- g0$par

draw_dgp <- function(n, seed) {
  gen_surv_weibull(n = n, k = k, pz = pz, q = q, SigU = par$SigU, A = par$A,
           beta = par$beta, alpha = par$alpha, B = par$B, sZ = par$sZ,
           sW = par$sW, T_par = par$T, C_par = par$C, T2_par = par$T2,
           pW = par$pW, seed = seed)
}

# inspect event/censoring proportions + pick evaluation times off the
# actual event-time distribution (avoid tail outliers as eval points)
g_check <- draw_dgp(1e5, 777777)
quants <- c(0.001, seq(0.01, 0.99, length.out = 99))
event_times <- g_check$data$event_time
fit_grid <- round(as.numeric(quantile(event_times, quants)), 10)
eval_grid <- round(as.numeric(quantile(event_times, c(0.25, 0.5, 0.75))), 10)


# checking interiority
print(prop.table(table(g_check$data$event)))
print(cif_interiority_check(g_check, eval_grid))

rm(g_check)
gc()

# independent population truth, evaluated over the full fit_grid so
# coverage can be plotted as a function of time, not just at eval_grid
g_truth <- draw_dgp(n_truth, 999999)
gt <- ground_truth_cr(g_truth, fit_grid, events)[effect %in% effects,
                                                 .(time_interest, event, effect, truth = value)]
gt_cells <- ground_truth_cr_cells(g_truth, fit_grid, events)
rm(g_truth)
gc()

# repeated estimation samples
#
# NOTE (A0/A1/A2 decomposition): requires modify_S = NULL (estimated CIF --
# otherwise there is no nuisance error to decompose), modify_G = "oracle"
# and martingale_debias = FALSE (so ria is the exact true-G IPCW core and
# the censoring term of E2 vanishes identically, making a1 a clean A1).
run_one <- function(seed) {
  data.table::setDTthreads(1)
  g <- draw_dgp(n, 10000 + seed)
  set.seed(20000 + seed)
  dr_obj <- one_step_debias_surv(as.data.table(g$data), X, Z, W, time_var,
                                 event_var, time_interest = fit_grid,
                                 martingale_debias = TRUE,
                                 # modify_G = "oracle",
                                 # modify_S = NULL,
                                 dgm = g)
  
  dr <- dr_obj$measures[event %in% events & effect %in% effects]
  dr[, `:=`(method = "DR", seed = seed, sample_size = n)]
  
  dr_po <- dr_obj$measures_po[event %in% events]
  dr_po[, `:=`(method = "DR", seed = seed, sample_size = n)]
  
  list(effects = dr, cells = dr_po)
}

run_out <- parallel::mclapply(seq_len(nseed), run_one, mc.cores = n_cores(),
                              mc.preschedule = FALSE, mc.set.seed = FALSE)
est_full <- rbindlist(lapply(run_out, `[[`, "effects"))
est_po_full <- rbindlist(lapply(run_out, `[[`, "cells"))

# pointwise operating characteristics
agg <- merge(est_full, gt, by = c("time_interest", "event", "effect"))
agg[, cov := value - 1.96 * sd <= truth & truth <= value + 1.96 * sd]

summary <- agg[, .(
  bias = mean(value - truth), empirical_sd = sd(value), mean_se = mean(sd),
  se_ratio = mean(sd) / sd(value), coverage = mean(cov),
  rmse = sqrt(mean((value - truth)^2))
), by = .(method, event, effect, time_interest)]

setorder(summary, event, effect, time_interest)
print(summary[time_interest %in% eval_grid])

# same pointwise operating characteristics, but for the raw psi(xz,xw,xy)
# cells (pre-differencing) -- helps isolate which nested counterfactual a
# bias comes from without re-deriving effects by hand
agg_po <- merge(est_po_full, gt_cells,
                by = c("xz", "xw", "xy", "event", "time_interest"))
agg_po[, cov := value - 1.96 * sd <= truth & truth <= value + 1.96 * sd]

summary_po <- agg_po[, .(
  bias = mean(value - truth), empirical_sd = sd(value), mean_se = mean(sd),
  se_ratio = mean(sd) / sd(value), coverage = mean(cov),
  rmse = sqrt(mean((value - truth)^2))
), by = .(method, event, xz, xw, xy, time_interest)]

setorder(summary_po, event, xz, xw, xy, time_interest)
print(summary_po[time_interest %in% eval_grid])

# ---- A0/A1/A2 decomposition (per seed, per cell) ----
# psi_hat - psi = A0 + A1 + A2 where
#   A0 = value_ora - truth    (empirical process, mean-zero Gaussian)
#   A2 = e2_t2 + e2_t3        (second-order remainder; censoring term = 0
#                              by oracle G, so this is the whole E2)
#   A1 = a1                   (cross-fitted empirical-process residual)
agg_po[, `:=`(
  A0 = value_ora - truth,
  A1 = a1,
  A2 = e2_t2 + e2_t3,
  bias_check = (value - truth) - (value_ora - truth) - a1 - (e2_t2 + e2_t3)
)]

# across-seed behavior on the sqrt(n)-scale, at the eval times
decomp <- agg_po[time_interest %in% eval_grid, .(
  bias = mean(value - truth),
  rn_A0_mean = mean(sqrt(sample_size) * A0),
  rn_A0_sd = sd(sqrt(sample_size) * A0),
  rn_A1_mean = mean(sqrt(sample_size) * A1),
  rn_A1_sd = sd(sqrt(sample_size) * A1),
  rn_E2t2_mean = mean(sqrt(sample_size) * e2_t2),
  rn_E2t2_sd = sd(sqrt(sample_size) * e2_t2),
  rn_E2t3_mean = mean(sqrt(sample_size) * e2_t3),
  rn_E2t3_sd = sd(sqrt(sample_size) * e2_t3),
  # A0 should be N(0, var(phi)): standardized by its own across-seed sd,
  # z0 quantiles should match N(0, 1)
  z0_p025 = quantile((A0 - mean(A0)) / sd(A0), 0.025),
  z0_p975 = quantile((A0 - mean(A0)) / sd(A0), 0.975),
  z0_mean_shift = mean(A0) / (sd(A0) / sqrt(.N))
), by = .(event, xz, xw, xy, time_interest)]

setorder(decomp, event, xz, xw, xy, time_interest)
print(decomp)

# headline read-out: which term carries the ctf-IE failure.
# ctfie = pso(0,0,1) - pso(0,1,1); expect e2_t2 dominant with opposite
# signs across the two cells
print(decomp[xz == 0 & xy == 1 & xw %in% c(0, 1) & event == 2])

saveRDS(list(estimates = est_full, truth = gt, summary = summary, agg = agg,
             estimates_po = est_po_full, truth_po = gt_cells, summary_po = summary_po,
             agg_po = agg_po, decomp = decomp, par = par),
        file = file.path(root, "results", f("cr-dml-weibull-{n}-zeta-96s.rds"))
        )

# bias_check should be ~0 by construction (up to NA-trimming); flag otherwise
stopifnot(max(abs(agg_po$bias_check), na.rm = TRUE) < 1e-8)