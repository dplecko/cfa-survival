
# pilot testbed: log(RMSE) vs log(n) across DR conditions
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

# ---- frozen DGM: same setup as synth-testbed.R ----
DGM_SEED <- 2026

pz <- 5
q <- 3
k <- 4

set.seed(DGM_SEED)
g0 <- gen_surv(n = 1, k = k, pz = pz, q = q, seed = DGM_SEED)
par <- g0$par

draw_dgp <- function(n, seed) {
  gen_surv(
    n = n, k = k, pz = pz, q = q,
    SigU = par$SigU, A = par$A, beta = par$beta,
    alpha = par$alpha, B = par$B,
    sZ = par$sZ, sW = par$sW,
    T_par = par$T, C_par = par$C,
    pW = par$pW, seed = seed
  )
}

effects <- c("tv", "ctfde", "ctfie", "ctfse")

# choose evaluation times and integration grid, quantile-based off the
# actual event-time distribution (avoids uniform-in-time tail outliers)
g_ref <- draw_dgp(n = 1e5, seed = 7777)
event_times <- g_ref$data$event_time
quants <- c(0.001, seq(0.01, 0.99, by = 0.01)) # 100 pts; hits 0.25/0.5/0.75 exactly

q_times <- round(as.numeric(quantile(event_times, c(0.25, 0.5, 0.75))), 10)
fine_grid <- round(as.numeric(quantile(event_times, quants)), 10)

rm(g_ref)
gc()

# independent population truth
g_truth <- draw_dgp(n = 1e6, seed = 888888)
gt <- ground_truth(g_truth, q_times)[
  effect %in% effects,
  .(effect, time_interest, gt = value)
]
rm(g_truth)
gc()

# ---- DR conditions ----
conditions <- list(
  list(name = "estimated", modify_S = NULL, modify_G = NULL),
  list(name = "S_wrong", modify_S = "corrupt", modify_G = "oracle"),
  list(name = "G_wrong", modify_S = "oracle", modify_G = "corrupt"),
  list(name = "both_wrong", modify_S = "corrupt", modify_G = "corrupt")
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
  
  g  <- draw_dgp(n = n_i, seed = sd_i)
  df <- as.data.table(g$data)
  
  learner_seed <- 1e6L + 1000L * match(n_i, n_grid) + sd_i
  
  rows <- list()
  for (cnd in conditions) {
    set.seed(learner_seed)
    dr_obj <- tryCatch(
      one_step_debias_surv(
        df, X, Z, W, time_var = "event_time", event_var = "event",
        time_interest = fine_grid,
        martingale_debias = TRUE,
        modify_S = cnd$modify_S,
        modify_G = cnd$modify_G,
        dgm = g
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

# analysis
conds <- c("estimated", "S_wrong", "G_wrong")

dr <- merge(
  res[condition %chin% conds & effect != "tv"],
  gt,
  by = c("effect", "time_interest")
)

diag <- dr[, .(
  bias = mean(value - gt),
  empirical_sd = sd(value),
  rmse = sqrt(mean((value - gt)^2))
), by = .(condition, effect, time_interest, n)]

diag[, `:=`(
  condition = factor(condition, levels = conds),
  t_lbl = sprintf("t = %.2f", time_interest)
)]

slopes_all <- diag[, .(
  slope_all = unname(coef(lm(log(rmse) ~ log(n)))[2])
), by = .(condition, effect, time_interest)]

slopes_tail <- diag[n >= 2000, .(
  slope_tail = unname(coef(lm(log(rmse) ~ log(n)))[2])
), by = .(condition, effect, time_interest)]

slopes <- merge(slopes_all, slopes_tail,
                by = c("condition", "effect", "time_interest"))
setorder(slopes, condition, effect, time_interest)
print(slopes)


diag[, effect := factor(effect, levels = c("ctfde", "ctfie", "ctfse"), 
                        labels = c("Ctf-DE", "Ctf-IE", "Ctf-SE"))]

diag[, Setting := factor(condition, levels = c("estimated", "S_wrong", "G_wrong"), 
                        labels = c("S, G estimated", "S misspecified, G oracle", 
                                   "G misspecified, S oracle"))]

p_rmse <- ggplot(diag, aes(n, rmse, color = Setting, group = Setting)) +
  geom_line() + geom_point() +
  scale_x_log10(breaks = n_grid) + scale_y_log10() +
  facet_grid(effect ~ t_lbl, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)) +
  labs(x = "n", y = "RMSE")

p_bias <- ggplot(diag, aes(n, bias, color = Setting, group = Setting)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() + geom_point() +
  scale_x_log10(breaks = n_grid) +
  facet_grid(effect ~ t_lbl, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "n", y = "Bias")

print(p_rmse)
print(p_bias)

ggsave(file.path(out_path, "rmse.png"), p_rmse, width = 10, height = 6, dpi = 150)
# ggsave(file.path(out_path, "bias.png"), p_bias, width = 11, height = 8, dpi = 150)
# 
# saveRDS(list(diag = diag, slopes = slopes, gt = gt),
#         file.path(out_path, "dr_summary.rds"))

# conds <- c("mart_on_Swrong", "mart_on_Gwrong", "mart_on_clean")
# effs <- c("ctfde", "ctfie", "ctfse")
# 
# tab <- copy(diag[condition %chin% conds & effect %chin% effs])
# tab[, `:=`(
#   condition = factor(condition, levels = conds),
#   effect = factor(effect, levels = effs,
#                   labels = c("Ctf-DE", "Ctf-IE", "Ctf-SE")),
#   entry = sprintf("%.2f", 100 * rmse)
# )]
# 
# cells <- tab[order(effect, time_interest, n, condition),
#              .(RMSE = paste(entry, collapse = " / ")),
#              by = .(effect, time_interest, n)]
# 
# cells[, row := sprintf("%s, $t=%.2f$", effect, time_interest)]
# out <- dcast(cells, row ~ n, value.var = "RMSE")
# 
# ns <- names(out)[-1]
# setnames(out, ns, paste0("$n=", ns, "$"))
# 
# knitr::kable(out, format = "pipe", align = c("l", rep("c", length(ns))))
