
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

effects <- c("tv", "ctfde", "ctfie", "ctfse")

# choose evaluation times and integration grid
g_ref <- draw(n = 1e5, seed = 7777)
q_times <- round(as.numeric(quantile(g_ref$data$T, c(0.25, 0.5, 0.75))), 10)
t_hi <- max(q_times, as.numeric(quantile(g_ref$data$event_time, 0.95)))
fine_grid <- sort(unique(round(c(seq(0.05, t_hi, length.out = 80), q_times), 10)))
rm(g_ref); gc()

# independent population truth
g_truth <- draw(n = 1e6, seed = 888888)
gt <- ground_truth(g_truth, q_times)[effect %in% effects, .(effect, time_interest, gt = value)]
rm(g_truth); gc()

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
  
  learner_seed <- 1e6L + 1000L * match(n_i, n_grid) + sd_i
  
  rows <- list()
  for (cnd in conditions) {
    set.seed(learner_seed)
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
    
    est <- dr_obj$measures[effect %in% effects]
    est[, time_interest := round(time_interest, 10)]
    est_q <- est[time_interest %in% q_times]
    stopifnot(nrow(est_q) == length(effects) * length(q_times))
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

invisible(mclapply(seq_len(nrow(work)), run_one, mc.cores = NCORES,
                   mc.preschedule = FALSE, mc.set.seed = FALSE))

part_files <- list.files(out_path, pattern = "^part_n.*\\.rds$", full.names = TRUE)
stopifnot(length(part_files) == nrow(work))
res <- rbindlist(lapply(part_files, readRDS))
stopifnot(nrow(res) == nrow(work) * length(conditions) * length(effects) * length(q_times))

# ---- plot ----
# obj <- readRDS(file.path(out_path, "pilot_summary.rds"))
# rmse <- obj$rmse
# rmse[, t_lbl := factor(sprintf("t=%.2f", time_interest))]
# p <- ggplot(rmse, aes(x = log_n, y = log_rmse, color = condition)) +
#   geom_line() + geom_point() +
#   facet_grid(effect ~ t_lbl, scales = "free_y") +
#   theme_bw() +
#   labs(x = "log n", y = "log RMSE",
#        title = "DR verification: log-RMSE vs log-n",
#        subtitle = sprintf("frozen DGM (seed %d), %d seeds/n", DGM_SEED, n_seeds))
# ggsave(file.path(out_path, "pilot_logrmse.png"),
#        plot = p, width = 11, height = 8)
# 
# print(obj$slopes)

# library(ggplot2); library(data.table)
# 
# obj  <- readRDS(file.path(out_path, "pilot_summary.rds"))
# rmse <- copy(obj$rmse)
# 
# rmse[, augmentation := fifelse(grepl("^mart_on", condition),
#                                "DR (augmented)", "IPCW only")]
# rmse[, nuisance := fcase(
#   grepl("clean",     condition), "S ok, G ok",
#   grepl("Swrong",    condition), "S wrong, G ok",
#   grepl("Gwrong",    condition) & !grepl("both", condition), "S ok, G wrong",
#   grepl("bothwrong", condition), "S wrong, G wrong"
# )]
# rmse[, nuisance := factor(nuisance, levels = c(
#   "S ok, G ok","S wrong, G ok","S ok, G wrong","S wrong, G wrong"))]
# rmse[, augmentation := factor(augmentation,
#                               levels = c("IPCW only","DR (augmented)"))]
# rmse[, effect := factor(effect, levels = c("tv","ctfde","ctfie","ctfse"),
#                         labels = c("TV","Ctf-DE","Ctf-IE","Ctf-SE"))]
# 
# # drop the late-time quantile — IPCW variance explodes; story is told at t1, t2
# # rmse <- rmse[time_interest <= quantile(rmse$time_interest, 0.67)]
# rmse[, t_lbl := sprintf("t = %.2f", time_interest)]
# 
# p <- ggplot(rmse, aes(x = log_n, y = log_rmse,
#                       color = nuisance, linetype = augmentation,
#                       shape = augmentation,
#                       group = interaction(nuisance, augmentation))) +
#   geom_line(linewidth = 0.7) +
#   geom_point(size = 2) +
#   facet_wrap(effect ~ t_lbl, scales = "free_y", ncol = 3) +
#   scale_color_manual("Nuisance status",
#                      values = c("S ok, G ok"       = "#1b9e77",   # green
#                                 "S wrong, G ok"    = "#e6ab02",   # darker yellow/gold
#                                 "S ok, G wrong"    = "#f4d03f",   # lighter yellow
#                                 "S wrong, G wrong" = "#d62728")) +
#   scale_linetype_manual("Augmentation", values = c("dashed","solid")) +
#   scale_shape_manual("Augmentation", values = c(1, 16)) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "bottom", legend.box = "vertical",
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 9)) +
#   labs(x = expression(log(n)), y = expression(log(RMSE)),
#        title = "Double-robustness verification",
#        subtitle = "log-RMSE vs log-n, frozen DGM, 64 seeds per n")
# 
# ggsave(file.path(out_path, "pilot_logrmse.png"),
#        plot = p, width = 11, height = 12, dpi = 150)
