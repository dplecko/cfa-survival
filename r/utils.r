
# file name helpers
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
  ggplot(km_dt[time < max_time], aes(x = time, y = surv, color = Majority, 
                                     fill = Majority)) +
    geom_line(linewidth=1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, linewidth=0) +
    theme_bw() + scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7),
          legend.box.background = element_rect()) + xlab("Time (days)") +
    xlim(c(0, 545))
}

tv_comparison <- function(fsurv, data) {
  
  sfit <- survfit(Surv(event_time, event) ~ majority, data = data)
  tgrid <- sort(unique(data[majority == 1]$event_time))
  sfit_sum <- summary(sfit, times = tgrid)
  
  s0 <- sfit_sum$surv[sfit_sum$strata == "majority=0"]
  s1 <- sfit_sum$surv[sfit_sum$strata == "majority=1"]
  
  km_dt <- data.table(
    surv = sfit_sum$surv, lwr = sfit_sum$lower, upr = sfit_sum$upper, 
    time = sfit_sum$time,
    Majority = as.integer(sfit_sum$strata) - 1
  )
  
  # determine the time after which KM estimate becomes unreliable
  max_time <- min(sfit$time[sfit$n.risk < round(0.005 * nrow(data))])
  
  if (!is.element("sd", fsurv$measures)) {
    
    fsurv$measures[, sd := 0]
  }
  
  ggplot(
    fsurv$measures[effect == "tv" & scale == "surv" & time_interest < max_time],
    aes(x = time_interest, y = value)
  ) +
    geom_line(linewidth=1) +
    geom_ribbon(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd),
                alpha = 0.2, linewidth=0) +
    theme_bw() + scale_y_continuous(labels = scales::percent) +
    # add TV from KM
    # geom_line(mapping = aes(x = time, y = V1),
    #           data = km_dt[time < max_time, diff(surv), by = c("time")],
    #           linewidth=1, color = "red", alpha = 0.6) +
    ylab("TV measure for survival") + xlab("Time (days)")
}

# compare TV decomposition with a static decomposition at T = t
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

# helper for survival curves and TV over time
compare_srv <- function(data, fsurv, type = c("tv", "marginal")) {
  
  type <- match.arg(type, c("tv", "marginal"))
  
  if (type == "marginal") {
    
    sfit <- survfit(Surv(event_time, event) ~ 1, data = data)
    ridx <- sfit$n.risk > 500
    
    dt_km <- data.table(surv = sfit$surv[ridx], time = sfit$time[ridx], 
                        method = "KM")
    
    dt_kmob <- km_oob(data)
    
    dt_rf <- data.table(surv = colMeans(fsurv$srv$srvx), 
                        time = fsurv$time_interest, method = "RF")
    
    
    p <- ggplot(rbind(dt_km, dt_kmob, dt_rf), 
                aes(x = time, y = surv, color = method)) +
      geom_line(linewidth=1) + 
      theme_bw() + scale_y_continuous(labels = scales::percent) +
      theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7),
            legend.box.background = element_rect()) + 
      xlab("Time (days)")
  } else if (type == "tv") {
    
    sfit <- survfit(Surv(event_time, event) ~ majority, data = data)
    dt_km <- data.table(
      surv = sfit$surv, time = sfit$time,
      majority = c(rep(0, sfit$strata[1]), rep(1, sfit$strata[2]))
    )
    # 
    dt_rf <- rbind(
      data.table(surv = colMeans(fsurv$srv$srvx[dat$majority == 0, ]),
                 time = fsurv$time_interest, majority = 0),
      data.table(surv = colMeans(fsurv$srv$srvx[dat$majority == 1, ]),
                 time = fsurv$time_interest, majority = 1)
    )
    
    dt <- rbind(dt_km[, method := "KM"], dt_rf[, method := "RF"])
    p1 <- ggplot(dt,
                 aes(x=time, y=surv, color = factor(majority), 
                     linetype = factor(method))) +
      geom_line(linewidth=1) + theme_bw() +
      coord_cartesian(xlim=c(0, 500), ylim = c(0.89, 0.95)) +
      scale_color_discrete(name = "Majority") +
      scale_linetype_manual(name = "Inference Method",
                            values = c("solid", "dotted")) +
      theme(legend.position = "inside", legend.position.inside = c(0.65, 0.65),
            legend.box.background = element_rect())
    
    p2 <- ggplot(
      data = dt[, diff(surv), by = c("time", "method")],
      aes(x = time, y = V1, color = method)
    ) +
      geom_line(linewidth=1) + 
      theme_bw() + scale_y_continuous(labels = scales::percent) +
      # add TV from KM
      geom_line(mapping = aes(x = time, y = V1), 
                data = dt_km[, diff(surv), by = c("time")],
                linewidth=0.8, color = "red", alpha = 0.6) +
      ylab("TV measure for survival") + xlab("Time (days)")
    
    p <- cowplot::plot_grid(p1, p2, ncol = 2L)
  }
  
  p
}
