
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# specify model parameters
pz <- 5
q <- 3
n <- 10^3
nseed <- 10

# construct the SFM
X <- "majority"
Z <- paste0("z", seq_len(pz))
W <- paste0("w", seq_len(q))
event_var <- "event"
time_var <- "event_time"

est_full <- gt_full <- c()
df_gen <- list()
for (seed in seq_len(nseed)) {
  
  # generate the data
  g <- gen_surv(n=n, k=4, pz=pz, q=q, seed=seed)
  df <- g$data
  df_gen[[seed]] <- df
  par <- g$par
  
  # get model-based estimates
  fsurv <- fair_surv(
    as.data.table(df), X, Z, W, time_var, event_var, time_interest = 100,
    nboot = 10
  )
  mod_est <- fsurv$measures[effect %in% c("tv", "ctfde", "ctfie", "ctfse") &
                            scale == "surv"]
  mod_est[, `:=`(method = "model", seed = seed, sample_size = n)]
  tgrid <- fsurv$time_interest
  
  # get DR estimates
  dr_obj <- one_step_debias_surv(as.data.table(df), X, Z, W, time_var, event_var, tgrid)
  dr_est <- dr_obj$measures
  dr_est[, `:=`(method = "DR", seed = seed, sample_size = n, scale = "surv",
                event = 1)]
  
  est_full <- rbind(est_full, mod_est, dr_est)
  
  # get ground truth
  gt_meas <- ground_truth(g, tgrid)[effect %in% c("tv", "ctfde", "ctfie", "ctfse")]
  gt_meas[, `:=`(method = "Truth", seed = seed)]
  gt_full <- rbind(gt_full, gt_meas)
}

# visual inspection for a single generative model
est_full[, method := factor(method)]
for (seed_num in seq_len(seed)) {
  
  p_curr <- ggplot(est_full[seed == seed_num & time_interest < 5], 
                   aes(x = time_interest, y = value, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd, fill = method),
                alpha = 0.4) +
    geom_line(data = gt_full[seed == seed_num], 
              color = "black") +
    coord_cartesian(xlim = c(0, quantile(df_gen[[seed]]$event_time, 0.9))) +
    facet_wrap(~ effect, scales = "free") + theme_bw()
  ggsave(plot = p_curr, 
         filename = file.path("results", "synth-tests", paste0("single-dgm-", seed_num, ".png")),
         width = 6, height = 6)
}

# aggregate coverage
by_vars <- c("time_interest", "effect", "seed")
agg_cov <- merge(est_full, gt_full[, c("value", by_vars), with=F], by = by_vars)
agg_cov <- agg_cov[time_interest < 5]
agg_cov[, cov := ((value.x + 1.96 * sd) > value.y) & 
                 ((value.x - 1.96 * sd) < value.y)]

ggplot(
  agg_cov[, list(coverage = mean(cov)), by = c("method", "effect")],
  aes(x = effect, y = coverage, fill = effect)
) +
  geom_col() + theme_bw() +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth=1) +
  facet_wrap(~ method)

