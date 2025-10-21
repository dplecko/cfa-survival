
fname <- function(src, out, balance, split) {
  
  balance <- if (!balance) NULL else "balanced"
  split <- if (!split) NULL else "split"
  paste0(paste(c(src, "fsurv", out, balance, split), collapse = "_"), ".RData")
}

km_curves <- function(data) {
  
  sfit <- survfit(Surv(event_time, event) ~ majority, data = data)
  
  # determine the time after which KM estimate becomes unreliable
  max_time <- min(sfit$time[sfit$n.risk < 100])
  
  km_dt <- data.table(
    surv = sfit$surv, lwr = sfit$lower, upr = sfit$upper, time = sfit$time,
    Majority = factor(c(rep(0, sfit$strata[1]), rep(1, sfit$strata[2])))
  )
  ggplot(km_dt[time < max_time], aes(x = time, y = surv, color = Majority, fill = Majority)) +
    geom_line(linewidth=1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, linewidth=0) +
    theme_bw() + scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7),
          legend.box.background = element_rect()) + xlab("Time (days)") +
    xlim(c(0, 545))
}

tv_comparison <- function(fsurv, data) {
  
  sfit <- survfit(Surv(event_time, event) ~ majority, data = data)
  
  # determine the time after which KM estimate becomes unreliable
  max_time <- min(sfit$time[sfit$n.risk < 100])
  
  km_dt <- data.table(
    surv = sfit$surv, lwr = sfit$lower, upr = sfit$upper, time = sfit$time,
    Majority = factor(c(rep(0, sfit$strata[1]), rep(1, sfit$strata[2])))
  )
  
  ggplot(
    fsurv$measures[effect == "tv" & scale == "surv"],
    aes(x = time_interest, y = value)
  ) +
    geom_line(linewidth=1) +
    geom_ribbon(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd),
                alpha = 0.2, linewidth=0) +
    theme_bw() + scale_y_continuous(labels = scales::percent) +
    # add TV from KM
    geom_line(mapping = aes(x = time, y = V1),
              data = km_dt[time < max_time, diff(surv), by = c("time")],
              linewidth=0.8, color = "red", alpha = 0.6) +
    ylab("TV measure for survival") + xlab("Time (days)")
}


srv_with_stat <- function(fsurv, data) {
  
  
  d_max <- as.Date("2024-07-01")
  T_hor <- 180
  d_max - T_hor
  
  dates <- load_concepts("icu_adm_date", "anzics")
  dates[, icu_adm_date := as.Date(icu_adm_date)]
  data <- merge(data, dates, all.x = TRUE)
  
  # filter!
  dat_stat <- data[icu_adm_date <= d_max - T_hor]
  dat_stat[event_time > T_hor & event == 1, y := 1]
  dat_stat[event_time <= T_hor & event == 1, y := 0]
  dat_stat[event == 0, y := 1]
  fcb <- faircause::fairness_cookbook(dat_stat, X = X, Z = Z, W = W, Y = "y", 
                                      x0 = 0, x1 = 1)
  
  # load("data/fcb_stat.RData")
  
  stat_plt <- cbind(fcb$measures, time = T_hor)
  stat_plt <- stat_plt[stat_plt$measure %in% c("ctfde", "ctfie", "ctfse"), ]
  names(stat_plt)[1] <- "effect"
  
  autoplot(fsurv, scale = "surv") +
    geom_point(data = stat_plt, aes(x = time, y = value)) +
    geom_errorbar(data = stat_plt, aes(x = time, y =  value,
                                       ymax = value + 1.96 * sd,
                                       ymin = value - 1.96 * sd), width=20) +
    # coord_cartesian(xlim = c(0, 200), ylim = c(-0.05, 0.025)) +
    scale_y_continuous(labels = scales::percent) + xlab("Time (days)") +
    theme(legend.position = "inside", legend.position.inside = c(0.2, 0.2),
          legend.box.background = element_rect())
}