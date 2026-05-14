
# pilot testbed: log(RMSE) vs log(n) across DR conditions
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

# ---- frozen DGM (parameters drawn once, reused for every sample) ----
DGM_SEED <- 42
pz <- 5; q <- 3; k <- 4

set.seed(DGM_SEED)
dgm <- list(
  k = k, pz = pz, q = q,
  SigU = diag(k),
  A     = matrix(rnorm(pz*k), pz, k),
  beta  = rnorm(k),
  alpha = rep(0.5, q),
  B     = matrix(rnorm(pz*q), pz, q),
  sZ    = rep(0.7, pz),
  sW    = rep(0.7, q),
  T_par = list(mu0=0,   muX= 0.6, muZ=rnorm(pz,0,0.2),
               muW=rnorm(q,0.2,0.2), sT=0.6),
  C_par = list(mu0=0.3, muX=-0.2, muZ=rnorm(pz,0,0.2),
               muW=rnorm(q,0.2,0.2), sT=1, winsL=0),
  pW    = runif(1, 0.6, 0.75)  # frozen, no longer drift across n
)

draw <- function(n, seed) {
  gen_surv(n=n, k=dgm$k, pz=dgm$pz, q=dgm$q,
           SigU=dgm$SigU, A=dgm$A, beta=dgm$beta,
           alpha=dgm$alpha, B=dgm$B,
           sZ=dgm$sZ, sW=dgm$sW,
           T_par=dgm$T_par, C_par=dgm$C_par,
           pW=dgm$pW, seed=seed)
}

# ---- eval times: 25/50/75% of marginal T (from large reference) ----
g_ref  <- draw(n = 1e5, seed = 7777)
q_times <- as.numeric(quantile(g_ref$data$T, c(0.25, 0.5, 0.75)))

# fine grid for the xi_2 Riemann integral; include q_times exactly
t_hi      <- as.numeric(quantile(g_ref$data$event_time, 0.95))
fine_grid <- sort(unique(c(seq(0.05, t_hi, length.out = 80), q_times)))

# ---- population ground truth (computed once on the reference) ----
gt <- ground_truth(g_ref, q_times)
gt <- gt[effect %in% c("tv","ctfde","ctfie","ctfse"),
         .(effect, time_interest, gt = value)]

# ---- conditions ----
conditions <- list(
  list(name="mart_off_clean",    mart=FALSE, cS=FALSE, cG=FALSE),
  list(name="mart_off_Gwrong",   mart=FALSE, cS=FALSE, cG=TRUE),
  list(name="mart_on_clean",     mart=TRUE,  cS=FALSE, cG=FALSE),
  list(name="mart_on_Swrong",    mart=TRUE,  cS=TRUE,  cG=FALSE),
  list(name="mart_on_Gwrong",    mart=TRUE,  cS=FALSE, cG=TRUE),
  list(name="mart_on_bothwrong", mart=TRUE,  cS=TRUE,  cG=TRUE)
)

# ---- pilot sweep ----
n_grid  <- c(500, 1000, 2000, 4000, 8000) # 
n_seeds <- 64 # 20

X <- "majority"
Z <- paste0("z", seq_len(pz))
W <- paste0("w", seq_len(q))

out_path <- file.path(root, "results", "dr-verification")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

results <- list()
library(parallel)

# detect cores from scheduler; UGE on Hoffman2 sets NSLOTS
NCORES <- n_cores()
cat("Using", NCORES, "cores\n")

# flat work list
work <- as.data.table(expand.grid(n = n_grid, seed = seq_len(n_seeds),
                                  KEEP.OUT.ATTRS = FALSE))

run_one <- function(idx) {
  # cap inner threading inside the worker
  data.table::setDTthreads(1)
  
  n_i  <- work$n[idx]
  sd_i <- work$seed[idx]
  
  task_file <- file.path(out_path, sprintf("part_n%05d_seed%03d.rds", n_i, sd_i))
  if (file.exists(task_file)) return(readRDS(task_file))  # resume support
  
  g  <- draw(n = n_i, seed = sd_i)
  df <- as.data.table(g$data)
  
  rows <- list()
  for (cnd in conditions) {
    dr_obj <- tryCatch(
      one_step_debias_surv(
        df, X, Z, W, time_var = "event_time", event_var = "event",
        time_interest = fine_grid,
        martingale_debias = cnd$mart,
        corrupt_S = cnd$cS, corrupt_G = cnd$cG
      ),
      error = function(e) { warning(conditionMessage(e)); NULL }
    )
    if (is.null(dr_obj)) next
    
    est <- dr_obj$measures[effect %in% c("tv","ctfde","ctfie","ctfse")]
    est_q <- rbindlist(lapply(q_times, function(tt) {
      sub <- est[, .SD[which.min(abs(time_interest - tt))], by = effect]
      sub[, time_interest := tt]; sub
    }))
    est_q[, `:=`(n = n_i, seed = sd_i, condition = cnd$name)]
    rows[[length(rows) + 1]] <-
      est_q[, .(n, seed, condition, effect, time_interest, value)]
  }
  
  bind <- rbindlist(rows)
  saveRDS(bind, task_file)
  cat(sprintf("[%s] done n=%d seed=%d\n",
              format(Sys.time(),"%H:%M:%S"), n_i, sd_i))
  bind
}

res_list <- mclapply(
  seq_len(nrow(work)), run_one,
  mc.cores       = NCORES,
  mc.preschedule = FALSE,   # uneven task durations → schedule on demand
  mc.set.seed    = FALSE    # we set seed via draw(), don't want mc to override
)
res <- rbindlist(res_list)

# ---- used if the chunks are already available ----
# res <- rbindlist(lapply(
#   list.files(out_path, pattern = "^part_n.*\\.rds$", full.names = TRUE),
#   readRDS
# ))

# ---- RMSE and slopes ----
res  <- merge(res, gt, by = c("effect","time_interest"))
res[, sq_err := (value - gt)^2]

rmse <- res[, .(rmse = sqrt(mean(sq_err)), nrep = .N),
            by = .(n, condition, effect, time_interest)]
rmse[, `:=`(log_n = log(n), log_rmse = log(rmse))]

slopes <- rmse[, .(slope = coef(lm(log_rmse ~ log_n))[2]),
               by = .(condition, effect, time_interest)]

saveRDS(list(raw = res, rmse = rmse, slopes = slopes, gt = gt,
             q_times = q_times, dgm = dgm),
        file.path(out_path, "pilot_summary.rds"))

# ---- plot ----
rmse[, t_lbl := factor(sprintf("t=%.2f", time_interest))]
p <- ggplot(rmse, aes(x = log_n, y = log_rmse, color = condition)) +
  geom_line() + geom_point() +
  facet_grid(effect ~ t_lbl, scales = "free_y") +
  theme_bw() +
  labs(x = "log n", y = "log RMSE",
       title = "DR verification: log-RMSE vs log-n",
       subtitle = sprintf("frozen DGM (seed %d), %d seeds/n", DGM_SEED, n_seeds))
ggsave(file.path(out_path, "pilot_logrmse.png"),
       plot = p, width = 11, height = 8)

print(slopes)