root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

# Informative censoring with known truth: competing-risks Weibull DGM (as in
# scripts/cr-synth.R) whose latent event times T1, T2 are coupled by a Frank
# copula given (X, Z, W), with Kendall's tau in {0.1, 0.5, 0.8}. Effects on the
# latent T1 are estimated with copula-graphic Routes I (stabilized, r1_cap = 10)
# and II, using the Frank copula at the true tau; n = 5000, 100 reps per tau.
#
#   SIM_LIST=1   Rscript scripts/ic-synth.R   # print the task table
#   SIM_TASK=<k> Rscript scripts/ic-synth.R   # run task k (or Slurm array id)
#   Rscript scripts/ic-synth.R agg            # truth + summary table
#   SIM_SMOKE=1  ...                          # smoke: 1 rep per tau, n = 1000

smoke <- nzchar(Sys.getenv("SIM_SMOKE"))
base <- file.path(root, "results", "ic-synth")
out_dir <- file.path(base, if (smoke) "smoke" else "parts")
X <- "majority"; ev <- "event"; tm <- "event_time"
DGM_SEED <- 2026

# ---------------------------------------------------------------- DGM -----
cr_dgm <- local({
  pz <- 2; q <- 2; k <- 2
  set.seed(DGM_SEED)
  T_par  <- list(mu0 = -0.15, muX = -0.15, muZ = c(-0.08, 0.07),
                 muW = c(-0.11, -0.15), shape = 1.5, scale = 1)
  T2_par <- list(mu0 = -0.15, muX =  0.15, muZ = c( 0.07, -0.06),
                 muW = c( 0.13,  0.13), shape = 1.5, scale = 1)
  C_par  <- list(mu0 = -0.85, muX = -0.08, muZ = c(-0.04, 0.04),
                 muW = c(-0.03, 0.03), shape = 1.2, scale = 1)
  par <- gen_surv_weibull(n = 1, k = k, pz = pz, q = q, C_par = C_par,
                          T_par = T_par, T2_par = T2_par, seed = DGM_SEED,
                          pW = 0.75, alpha = rep(0.75, q))$par
  list(Z = paste0("z", seq_len(pz)), W = paste0("w", seq_len(q)),
       draw = function(n, seed, copula_tau = NULL)
         gen_surv_weibull(n = n, k = k, pz = pz, q = q, SigU = par$SigU,
                          A = par$A, beta = par$beta, alpha = par$alpha,
                          B = par$B, sZ = par$sZ, sW = par$sW, T_par = par$T,
                          C_par = par$C, T2_par = par$T2, pW = par$pW,
                          seed = seed, copula_tau = copula_tau))
})

# 12 equally spaced points up to the 90% quantile of M (tau = 0.5 draw);
# evaluated at grid points 3, 6, 9 (~25/50/75% of the grid)
time_grid <- function() {
  h <- as.numeric(quantile(cr_dgm$draw(1e5, 777777, copula_tau = 0.5)$data$event_time, 0.9))
  tg <- round(h * seq_len(12) / 12, 10)
  list(fit = tg, eval = tg[c(3, 6, 9)])
}

taus <- c(0.1, 0.5, 0.8)
n_obs <- if (smoke) 1000 else 5000
n_rep <- if (smoke) 1 else 100
blk <- as.integer(Sys.getenv("SIM_BLK", "8"))   # reps per task
tasks <- do.call(rbind, lapply(taus, function(tau) {
  st <- seq(1, n_rep, by = blk)
  data.frame(tau = tau, from = st, to = pmin(st + blk - 1, n_rep))
}))
tasks$task <- seq_len(nrow(tasks))

# ---------------------------------------------------------------- agg -----
if (identical(commandArgs(TRUE)[1], "agg")) {
  # truth: counterfactual survival of the latent T1 per world (closed form, no
  # copula needed), averaged over 8 independent draws of 1e6
  truth_file <- file.path(base, "truth.rds")
  if (!file.exists(truth_file)) {
    eff <- measure_spec(); gr <- time_grid()
    one <- function(s) {
      data.table::setDTthreads(1)
      cells <- ground_truth_cells(cr_dgm$draw(1e6, s), gr$eval)
      cl <- function(sp, t) cells[xz == sp[1] & xw == sp[2] & xy == sp[3] & time_interest == t, truth]
      rbindlist(lapply(eff, function(e) rbindlist(lapply(gr$eval, function(t)
        data.table(effect = e$nm, time_interest = t,
                   truth = cl(e$spc[[1]], t) - cl(e$spc[[2]], t))))))
    }
    tr <- rbindlist(parallel::mclapply(31000 + 1:8, one, mc.cores = 4,
                                       mc.preschedule = FALSE))
    tr <- tr[, .(truth = mean(truth), truth_mcse = sd(truth) / sqrt(.N)),
             by = .(effect, time_interest)]
    saveRDS(tr, truth_file)
  }
  tr <- readRDS(truth_file)[, time_interest := round(time_interest, 10)]

  res <- rbindlist(lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
  res[, time_interest := round(time_interest, 10)]
  res <- merge(res, tr, by = c("effect", "time_interest"))
  tg <- sort(unique(res$time_interest))
  res[, q := c("25%", "50%", "75%")[match(time_interest, tg)]]

  # unstable reps: identified set inverted or extreme SE at any effect/time
  cat("unstable reps (id_hi < id_lo or sd > 0.2) by route x tau:\n")
  print(dcast(res[, .(bad = any(id_hi < id_lo | sd > 0.2)), by = .(route, tau, rep)][
    , .(bad = sum(bad)), by = .(route, tau)], route ~ tau, value.var = "bad"))

  tab <- res[, .(reps = .N,
                 coverage = mean(lwr <= truth & truth <= upr),           # partial-ID CI
                 cov_midpointCI = mean(lwr_mid <= truth & truth <= upr_mid),
                 truth_in_idset = mean(id_lo <= truth & truth <= id_hi),
                 med_idset_w = median(id_hi - id_lo), med_pidCI_w = median(upr - lwr),
                 bias = mean(value - truth), rmse = sqrt(mean((value - truth)^2))),
             by = .(route, tau, effect, q)]
  setorder(tab, route, tau, effect, q)
  print(tab, digits = 3, nrows = 500)
  fwrite(tab, file.path(base, if (smoke) "table-smoke.csv" else "table.csv"))
  quit(save = "no")
}

if (nzchar(Sys.getenv("SIM_LIST"))) {
  print(tasks, row.names = FALSE)
  quit(save = "no")
}

# ---------------------------------------------------------------- run -----
k_env <- Sys.getenv("SIM_TASK", Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
stopifnot(nzchar(k_env))
tk <- tasks[as.integer(k_env), ]
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("SIM_CORES", "1")))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("task %d: tau %s reps %d-%d on %d cores\n", tk$task, tk$tau, tk$from, tk$to, cores))

one_rep <- function(r) {
  data.table::setDTthreads(1)
  t0 <- Sys.time()
  gr <- time_grid()
  g <- cr_dgm$draw(n_obs, 100000L + 1000L * match(tk$tau, taus) + r, copula_tau = tk$tau)
  set.seed(200000L + 1000L * match(tk$tau, taus) + r)
  # nuisances are estimated (dgm = NULL): the oracle code assumes independent
  # latent times, which does not hold under this DGM
  o <- one_step_debias_surv(as.data.table(g$data), X, cr_dgm$Z, cr_dgm$W, tm, ev,
                            time_interest = gr$fit, copula = "frank",
                            tau_grid = tk$tau, route = c("I", "II"))
  res <- o$measures_sens[time_interest %in% gr$eval]
  res[, `:=`(tau_true = tk$tau, n = n_obs, rep = r,
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))]
  res
}

out <- rbindlist(parallel::mclapply(
  tk$from:tk$to, function(r) tryCatch(one_rep(r), error = function(e) {
    message("rep ", r, " failed: ", conditionMessage(e)); NULL }),
  mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE), fill = TRUE)

f_out <- file.path(out_dir, sprintf("tau%s_r%03d-%03d.rds", tk$tau, tk$from, tk$to))
saveRDS(out, f_out)
cat("saved", f_out, "with", nrow(out), "rows; reps ok:", length(unique(out$rep)), "of",
    tk$to - tk$from + 1, "\n")
