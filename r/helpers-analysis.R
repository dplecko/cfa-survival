
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
