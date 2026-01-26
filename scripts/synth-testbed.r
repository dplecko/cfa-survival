
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# specify model parameters
pz <- 5
q <- 3
n <- 10^4
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
  dr_est <- one_step_debias_surv(as.data.table(df), X, Z, W, time_var, event_var, tgrid)
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

# ggplot(rbind(mod_est),
#        aes(x = time_interest, y = value, color = method)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd, fill = method),
#               alpha = 0.4) +
#   geom_line(data = gt_meas,
#             color = "black") +
#   coord_cartesian(xlim = c(0, quantile(df_gen[[seed]]$event_time, 0.9))) +
#   facet_wrap(~ effect, scales = "free") + theme_bw()

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

# fit fair surv models on the data
# fsurv_v2m <- fair_surv_v2(as.data.table(df), X, Z, W, time_var, event_var, 
#                           method = "model")
# 
# fsurv_v2i <- fair_surv_v2(as.data.table(df), X, Z, W, time_var, event_var, 
#                           method = "ipw")

# get the ground truth

# testing TV rules
gt_meas[effect %in% c("tv", "ctfde", "ctfie", "ctfse"), 
        value[1] - value[2] + value[3] + value[4], 
        by = c("time")]

ggplot(
  rbind(
    fsurv_v2m$measures[, method := "model-based"],
    fsurv_v2i$measures[, method := "ipw"],
    gt_meas[, method := "truth"]
  ), aes(x = time, y = value, color = method)
) + geom_line(linewidth=1) + theme_bw() +
  facet_wrap(~ effect) +
  xlim(c(0, quantile(df$event_time, probs = 0.95)))

# get xgb fit
xgb_obj <- xgb_surv_cf(df, X, time_var, event_var, rhs = ".", time_interest = 150,
                       balance_groups = FALSE, split_forest = FALSE)

plot(xgb_obj$time_interest, colMeans(xgb_obj$srv), pch = 19,
     ylim = c(0, 1))

sfit <- survfit(Surv(event_time, event) ~ 1, data = df)

lines(sfit$time, sfit$surv, col = "blue")


#' * fit the standard fair-surv -> is there a TV-faithfulness issue? *
fsurv_synth <- fair_surv(as.data.table(df), X, Z, W, time_var, event_var)

# tv_comparison(fsurv_synth, as.data.table(df))


# get a random forest fit with ranger
# rfit <- ranger(
#   formula = Surv(event_time, event) ~ majority + z1 + z2 + z3 + z4 + z5 + w1 + w2 + w3,
#   data = df
# )

# get the TV measure
# rangl <- reshape2::melt(rfit$survival)
# names(rangl) <- c("row", "time_idx", "surv")
# rangl <- as.data.table(rangl)
# 
# tt <- data.table(time = rfit$unique.death.times, time_idx = seq_along(rfit$unique.death.times))
# maj <- data.table(row = seq_len(nrow(df)), majority = df$majority)
# rangl <- merge(rangl, tt)
# rangl <- merge(rangl, maj, by = "row")
# rangr <- rangl[, list(surv = mean(surv)), by = c("time", "majority")]


tv_comparison(fsurv_synth, as.data.table(df)) +
# geom_line(mapping = aes(x = time, y = V1),
#           data = rangr[time < 6, diff(surv), by = c("time")],
#           linewidth=1, color = "blue", alpha = 0.6) +
  geom_line(
    mapping = aes(x = time, y = value),
    ground_truth(g, tg = seq(0, 6, length.out = 100))[effect == "tv"],
    color = "pink", linewidth = 1
  )
