
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

DGM_SEED <- 2026

# settings
pz <- 5
q <- 3
k <- 4
n <- 10^4
nseed <- 96
n_truth <- 1e6
effects <- c("tv", "ctfde", "ctfie", "ctfse")

X <- "majority"
Z <- paste0("z", seq_len(pz))
W <- paste0("w", seq_len(q))
event_var <- "event"
time_var <- "event_time"

# fix one DGP
set.seed(DGM_SEED)
g0 <- gen_surv(n = 1, k = k, pz = pz, q = q, seed = DGM_SEED)
par <- g0$par

draw_dgp <- function(n, seed) {
  gen_surv(n = n, k = k, pz = pz, q = q, SigU = par$SigU, A = par$A, beta = par$beta,
           alpha = par$alpha, B = par$B, sZ = par$sZ, sW = par$sW,
           T_par = par$T, C_par = par$C, pW = par$pW, seed = seed)
}

# quantile-based evaluation/fit grids off the actual event-time distribution
g_check <- draw_dgp(1e5, 777777)
quants <- c(0.001, seq(0.01, 0.99, by = 0.01)) # 100 pts; hits 0.25/0.5/0.75 exactly
event_times <- g_check$data$event_time
fit_grid  <- round(as.numeric(quantile(event_times, quants)), 10)
eval_grid <- round(as.numeric(quantile(event_times, c(0.25, 0.5, 0.75))), 10)
rm(g_check); gc()

# independent population truth
g_truth <- draw_dgp(n_truth, 999999)
gt <- ground_truth(g_truth, eval_grid)[effect %in% effects,
                                       .(time_interest, effect, truth = value)]
rm(g_truth)
gc()

# repeated estimation samples
est_full <- vector("list", nseed)

run_one <- function(seed) {
  data.table::setDTthreads(1)
  g <- draw_dgp(n, 10000 + seed)
  set.seed(20000 + seed)
  dr <- one_step_debias_surv(as.data.table(g$data), X, Z, W, time_var, event_var,
                             time_interest = fit_grid, dgm = g)$measures
  dr <- dr[effect %in% effects & time_interest %in% eval_grid]
  dr[, `:=`(method = "DR", seed = seed, sample_size = n, scale = "surv", event = 1)]
  dr
}

est_full <- rbindlist(
  parallel::mclapply(seq_len(nseed), run_one, mc.cores = n_cores(),
                     mc.preschedule = FALSE, mc.set.seed = FALSE)
)

# pointwise operating characteristics
agg <- merge(est_full, gt, by = c("time_interest", "effect"))
agg[, cov := value - 1.96 * sd <= truth & truth <= value + 1.96 * sd]

summary <- agg[, .(
  bias = mean(value - truth),
  empirical_sd = sd(value),
  mean_se = mean(sd),
  se_ratio = mean(sd) / sd(value),
  coverage = mean(cov)
), by = .(method, effect, time_interest)]

setorder(summary, effect, time_interest)
print(summary)

saveRDS(list(estimates = est_full, truth = gt, summary = summary, agg = agg, 
             par = par),
        file.path(root, "results", "nic-dml-coverage-10k.rds"))

s   <- readRDS(file.path(root, "results", "nic-dml-coverage-10k.rds"))$summary
tab <- dcast(s, effect ~ time_interest, value.var = "coverage")
for (j in names(tab)[-1]) tab[[j]] <- sprintf("%.0f\\%%", 100 * tab[[j]])

knitr::kable(
  tab, format = "latex", booktabs = TRUE, escape = FALSE,
  align = c("l", rep("r", ncol(tab) - 1L)),
  col.names = c("Effect", names(tab)[-1]),
  caption = "Coverage of DML confidence intervals."
)
