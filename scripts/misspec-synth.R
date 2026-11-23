root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

# Nuisance misspecification beyond scripts/nic-misspec.R. "oracle" = true
# nuisance, "corrupt" = marginal KM / Aalen-Johansen (S, G) or marginal P(x);
# the remaining nuisances (always nu) are estimated.
#   S2  (NIC DGM of scripts/nic-misspec.R): S oracle,  G and P(x) corrupt
#   C1a (CR DGM of scripts/cr-synth.R):     S corrupt, G oracle
#   C1b (CR):                               S oracle,  G corrupt
#   C2  (CR):                               S oracle,  G and P(x) corrupt
# n in {500, ..., 8000}, 64 seeds; Wald CIs, effects at 25/50/75% of M.
#
#   SIM_LIST=1   Rscript scripts/misspec-synth.R   # print the task table
#   SIM_TASK=<k> Rscript scripts/misspec-synth.R   # run task k (or Slurm array id)
#   Rscript scripts/misspec-synth.R agg            # truth + summary table
#   SIM_SMOKE=1  ...                               # smoke: 2 seeds, n = 500

smoke <- nzchar(Sys.getenv("SIM_SMOKE"))
base <- file.path(root, "results", "misspec-synth")
out_dir <- file.path(base, if (smoke) "smoke" else "parts")
X <- "majority"; ev <- "event"; tm <- "event_time"
DGM_SEED <- 2026

# ---------------------------------------------------------------- DGMs ----
nic_dgm <- local({
  pz <- 5; q <- 3; k <- 4
  set.seed(DGM_SEED)
  par <- gen_surv(n = 1, k = k, pz = pz, q = q, seed = DGM_SEED)$par
  list(Z = paste0("z", seq_len(pz)), W = paste0("w", seq_len(q)),
       draw = function(n, seed)
         gen_surv(n = n, k = k, pz = pz, q = q, SigU = par$SigU, A = par$A,
                  beta = par$beta, alpha = par$alpha, B = par$B, sZ = par$sZ,
                  sW = par$sW, T_par = par$T, C_par = par$C, pW = par$pW,
                  seed = seed))
})
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
       draw = function(n, seed)
         gen_surv_weibull(n = n, k = k, pz = pz, q = q, SigU = par$SigU,
                          A = par$A, beta = par$beta, alpha = par$alpha,
                          B = par$B, sZ = par$sZ, sW = par$sW, T_par = par$T,
                          C_par = par$C, T2_par = par$T2, pW = par$pW,
                          seed = seed))
})

# 100-point quantile fit grid of M, evaluated at its 25/50/75% quantiles
time_grid <- function(setting) {
  if (setting == "nic") {
    et <- nic_dgm$draw(1e5, 7777)$data$event_time
    quants <- c(0.001, seq(0.01, 0.99, by = 0.01))
  } else {
    et <- cr_dgm$draw(1e5, 777777)$data$event_time
    quants <- c(0.001, seq(0.01, 0.99, length.out = 99))
  }
  list(fit = round(as.numeric(quantile(et, quants)), 10),
       eval = round(as.numeric(quantile(et, c(0.25, 0.5, 0.75))), 10))
}

arms <- list(
  S2  = list(setting = "nic", modify_S = "oracle",  modify_G = "corrupt", modify_P = "corrupt"),
  C1a = list(setting = "cr",  modify_S = "corrupt", modify_G = "oracle",  modify_P = NULL),
  C1b = list(setting = "cr",  modify_S = "oracle",  modify_G = "corrupt", modify_P = NULL),
  C2  = list(setting = "cr",  modify_S = "oracle",  modify_G = "corrupt", modify_P = "corrupt")
)
n_grid <- if (smoke) 500 else c(500, 1000, 2000, 4000, 8000)
n_seeds <- if (smoke) 2 else 64
# seeds per task (roughly constant run time per task)
blk <- function(n) c(`500` = 64, `1000` = 32, `2000` = 16, `4000` = 8, `8000` = 8)[[as.character(n)]]
tasks <- do.call(rbind, lapply(names(arms), function(a) do.call(rbind, lapply(n_grid, function(n) {
  st <- seq(1, n_seeds, by = blk(n))
  data.frame(arm = a, n = n, from = st, to = pmin(st + blk(n) - 1, n_seeds))
}))))
tasks$task <- seq_len(nrow(tasks))

# ---------------------------------------------------------------- agg -----
if (identical(commandArgs(TRUE)[1], "agg")) {
  # truth averaged over 8 independent draws of 1e6 per DGM
  truth_file <- file.path(base, "truth.rds")
  if (!file.exists(truth_file)) {
    gN <- time_grid("nic"); gC <- time_grid("cr")
    one <- function(s) {
      data.table::setDTthreads(1)
      tC <- ground_truth_cr(cr_dgm$draw(1e6, s), gC$eval, 1:2)[
        , .(setting = "cr", effect, time_interest, event, truth = value)]
      gc()
      tN <- ground_truth(nic_dgm$draw(1e6, s + 100), gN$eval)[
        , .(setting = "nic", effect, time_interest, event = 1L, truth = value)]
      rbind(tC, tN)
    }
    tr <- rbindlist(parallel::mclapply(31000 + 1:8, one, mc.cores = 4,
                                       mc.preschedule = FALSE))
    tr <- tr[, .(truth = mean(truth), truth_mcse = sd(truth) / sqrt(.N)),
             by = .(setting, effect, time_interest, event)]
    saveRDS(tr, truth_file)
  }
  tr <- readRDS(truth_file)[, time_interest := round(time_interest, 10)]

  res <- rbindlist(lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
  res[, setting := ifelse(arm == "S2", "nic", "cr")]
  res <- merge(res, tr, by = c("setting", "effect", "time_interest", "event"))
  tab <- res[, .(reps = .N, bias = mean(value - truth), emp_sd = sd(value),
                 rmse = sqrt(mean((value - truth)^2)),
                 coverage = mean(abs(value - truth) <= 1.96 * sd),
                 truth_mcse = truth_mcse[1]),
             by = .(arm, event, effect, time_interest, n)]
  setorder(tab, arm, event, effect, time_interest, n)
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
cat(sprintf("task %d: arm %s n %d seeds %d-%d on %d cores\n",
            tk$task, tk$arm, tk$n, tk$from, tk$to, cores))

one_rep <- function(r) {
  data.table::setDTthreads(1)
  t0 <- Sys.time()
  a <- arms[[tk$arm]]; d <- if (a$setting == "nic") nic_dgm else cr_dgm
  gr <- time_grid(a$setting)
  g <- d$draw(tk$n, r)                                  # seeds as nic-misspec.R
  set.seed(1e6L + 1000L * match(tk$n, c(500, 1000, 2000, 4000, 8000)) + r)
  o <- one_step_debias_surv(as.data.table(g$data), X, d$Z, d$W, tm, ev,
                            time_interest = gr$fit, martingale_debias = TRUE,
                            modify_S = a$modify_S, modify_G = a$modify_G,
                            modify_P = a$modify_P, dgm = g)
  res <- o$measures[effect %in% c("tv", "ctfde", "ctfie", "ctfse")]
  res[, time_interest := round(time_interest, 10)]
  res <- res[time_interest %in% gr$eval]
  if (!"event" %in% names(res)) res[, event := 1L]
  res[, `:=`(arm = tk$arm, n = tk$n, rep = r,
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))]
  res
}

out <- rbindlist(parallel::mclapply(
  tk$from:tk$to, function(r) tryCatch(one_rep(r), error = function(e) {
    message("rep ", r, " failed: ", conditionMessage(e)); NULL }),
  mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE), fill = TRUE)

f_out <- file.path(out_dir, sprintf("%s_n%05d_r%03d-%03d.rds", tk$arm, tk$n, tk$from, tk$to))
saveRDS(out, f_out)
cat("saved", f_out, "with", nrow(out), "rows; reps ok:", length(unique(out$rep)), "of",
    tk$to - tk$from + 1, "\n")
