root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

DGM_SEED <- 2026

# settings
pz <- 5
q <- 3
k <- 4
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
# C_par_neutral <- list(mu0 = 0.55, muX = 0, muZ = rep(0, 5),
#   muW = rep(0, 3), sT = 1, winsL = 0)
T_par <- list(mu0 = 0.25, muX = 0.15,
  muZ = c(0.06, -0.05, 0.04, 0.03, -0.05),
  muW = c(0.05, -0.04, 0.05), sT = 1.10)

T2_par <- list(mu0 = 0.45, muX = -0.15,
  muZ = c(-0.05, 0.04, -0.03, 0.05, 0.03),
  muW = c(-0.04, 0.05, -0.04), sT = 1.10)

C_par <- list(mu0 = 0.65, muX = 0.08,
  muZ = c(0.03, -0.03, 0.02, 0.02, -0.02),
  muW = c(0.025, -0.02, 0.025), sT = 1.20, winsL = 0)
g0 <- gen_surv(n = 1, k = k, pz = pz, q = q, C_par = C_par,
  T_par = T_par, T2_par = T2_par, seed = DGM_SEED)
par <- g0$par

draw_dgp <- function(n, seed) {
  gen_surv(n = n, k = k, pz = pz, q = q, SigU = par$SigU, A = par$A,
    beta = par$beta, alpha = par$alpha, B = par$B, sZ = par$sZ,
    sW = par$sW, T_par = par$T, C_par = par$C, T2_par = par$T2,
    pW = par$pW, seed = seed)
}

# inspect event/censoring proportions + pick evaluation times off the
# actual event-time distribution (avoid tail outliers as eval points)
g_check <- draw_dgp(1e5, 777777)
print(prop.table(table(g_check$data$event)))
quants <- c(0.001, seq(0.01, 0.99, length.out = 99))
event_times <- g_check$data$event_time
fit_grid <- round(as.numeric(quantile(event_times, quants)), 10)
eval_grid <- round(as.numeric(quantile(event_times, c(0.25, 0.5, 0.75))), 10)

# cif_gt <- CIF_conditional_exact(g_check, eval_grid, event = 2)
# 
# ggplot(reshape2::melt(cif_gt$cifx0), aes(x = value)) +
#   stat_ecdf() + facet_wrap(~ Var2)
# 
# colMeans(cif_gt$cifx0 < 0.01)
# colMeans(cif_gt$cifx0 < 0.05)
# colMeans(cif_gt$cifx0 < 0.1)
# colMeans(cif_gt$cifx0 > 0.6)

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
run_one <- function(seed) {
  data.table::setDTthreads(1)
  g <- draw_dgp(n, 10000 + seed)
  set.seed(20000 + seed)
  dr_obj <- one_step_debias_surv(as.data.table(g$data), X, Z, W, time_var,
    event_var, time_interest = fit_grid, 
    martingale_debias = FALSE,
    modify_G = "oracle",
    modify_S = "oracle",
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

saveRDS(list(estimates = est_full, truth = gt, summary = summary, agg = agg,
  estimates_po = est_po_full, truth_po = gt_cells, summary_po = summary_po,
  agg_po = agg_po, par = par),
  file.path(root, "results", "cr-dml-coverage-zeta-nDGP-10k-fora-nozeta.rds"))

s_po <- readRDS(file.path(root, "results",
  "cr-dml-coverage-zeta-nDGP-10k-fora.rds"))$summary_po
tab <- dcast(s_po[time_interest %in% eval_grid],
  event + xz + xw + xy ~ time_interest, value.var = "coverage")
for (j in names(tab)[-(1:4)]) {
  tab[[j]] <- sprintf("%.0f%%", 100 * tab[[j]])
}
knitr::kable(tab, format = "pipe",
  align = c("r", "l", rep("r", length(eval_grid))))

s <- readRDS(file.path(root, "results",
                       "cr-dml-coverage-zeta-nDGP-10k-fora.rds"))$summary
tab <- dcast(s[time_interest %in% eval_grid],
             event + effect ~ time_interest, value.var = "coverage")
for (j in names(tab)[-(1:2)]) {
  tab[[j]] <- sprintf("%.0f%%", 100 * tab[[j]])
}
knitr::kable(tab, format = "pipe",
             align = c("r", "l", rep("r", length(eval_grid))))

# coverage as a function of time (full fit_grid, not just eval_grid)
p_cov <- ggplot(s, aes(time_interest, coverage, color = effect)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_vline(xintercept = eval_grid, linetype = "dotted", alpha = 0.4) +
  facet_wrap(~event) +
  theme_bw() +
  labs(x = "time", y = "coverage", title = "Coverage over time")
print(p_cov)
