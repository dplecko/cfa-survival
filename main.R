#!/burg/opt/R-4.3.1/bin/Rscript
#SBATCH --job-name=fairsurv
#SBATCH --account=dsi
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --time=24:00:00

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# select data source
src <- "aics"

# prepare the data and the SFM
dat <- load_data("aics")
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")

set.seed(2026)

# for testing
local <- FALSE
if (local) {
  
  dat_run <- rbind(
    dat[majority == 0],
    dat[sample(which(dat$majority == 1), size = sum(dat$majority == 0))]
  )
  # dat_run <- dat[sample.int(nrow(dat), size = 10000, replace=FALSE)]
  nboot <- 3
} else {
  
  dat_run <- dat
  nboot <- 64
}

fsurv <- fair_surv(dat_run, X, Z, W, time_var, event_var, 
                   method = "rfs-cf", nboot = nboot,
                   balance_groups = TRUE)

if (!local) save(fsurv, file = paste0("data/", src, "_fsurvb.RData"))

# load(paste0("data/", src, "_fsurvb.RData"))

if (local) {
  
  ### paper plots:
  
  # (A) Kaplan-Meier Survival Curve Estimates
  sfit <- survfit(Surv(event_time, event) ~ majority, data = dat_run)
  km_dt <- data.table(
    surv = sfit$surv, lwr = sfit$lower, upr = sfit$upper, time = sfit$time,
    Majority = factor(c(rep(0, sfit$strata[1]), rep(1, sfit$strata[2])))
  )
  ggplot(km_dt, aes(x = time, y = surv, color = Majority, fill = Majority)) +
    geom_line(linewidth=1) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, linewidth=0) + 
    theme_bw() + scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7),
          legend.box.background = element_rect()) + xlab("Time (days)")
  ggsave("results/surv-curves.png", width = 6, height = 4)
  
  # (B) Survival-TV estimate from fairsurv object
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
              data = km_dt[, diff(surv), by = c("time")],
              linewidth=0.8, color = "red", alpha = 0.6) +
    ylab("TV measure for survival") + xlab("Time (days)")
  ggsave("results/surv-tv-ot.png", width = 6, height = 4)
  
  # (C) Survival-TV decomposition with a static comparison
  d_max <- as.Date("2024-07-01")
  T_hor <- 180
  d_max - T_hor
  
  dates <- load_concepts("icu_adm_date", "anzics")[, 
                                                   icu_adm_date := as.Date(icu_adm_date)]
  dat_run <- merge(dat_run, dates, all.x = TRUE)
  
  # filter!
  dat_stat <- dat_run[icu_adm_date <= d_max - T_hor]
  dat_stat[event_time > T_hor & event == 1, y := 1]
  dat_stat[event_time <= T_hor & event == 1, y := 0]
  dat_stat[event == 0, y := 1]
  #
  library(faircause)
  fcb <- fairness_cookbook(dat_stat, X = X, Z = Z, W = W, Y = "y", x0 = 0, x1 = 1)
  
  # load("data/fcb_stat.RData")
  
  stat_plt <- cbind(fcb$measures, time = T_hor)
  stat_plt <- stat_plt[stat_plt$measure %in% c("ctfde", "ctfie", "ctfse"), ]
  names(stat_plt)[1] <- "effect"

  autoplot(fsurv, scale = "surv") +
    geom_point(data = stat_plt, aes(x = time, y = value)) +
    geom_errorbar(data = stat_plt, aes(x = time, y =  value,
                                       ymax = value + 1.96 * sd,
                                       ymin = value - 1.96 * sd), width=20) +
    coord_cartesian(xlim = c(0, 200), ylim = c(-0.05, 0.025)) + 
    scale_y_continuous(labels = scales::percent) + xlab("Time (days)") +
    theme(legend.position = "inside", legend.position.inside = c(0.2, 0.2),
          legend.box.background = element_rect())
  ggsave("results/decomp-with-static.png", width = 6, height = 4)
  
  # (D) CHF-TV decomposition with a Cox comparison
  autoplot(fsurv, scale = "chf-ratio") +
    coord_cartesian(ylim=c(0.75, 1.8)) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7),
          legend.box.background = element_rect())
  
  ggsave("results/chf-ratio.png", width = 6, height = 4)
}


#' # marginal survival curves: what happens? (underestimation?)
#' sfit <- survfit(Surv(event_time, event) ~ 1, data = dat)
#' 
#' plot(sfit, ylim = c(0.9, 1))
#' lines(fsurv$time_interest, colMeans(fsurv$srv$srvx), col = "orange")
#' #' * Kaplan-Meier vs. RSF marginal fit has a gap of about 0.16% *
#' 
#' # what happens group-wise?
#' sfit <- survfit(Surv(event_time, event) ~ majority, data = dat)
#' km_dt <- data.table(
#'   surv = sfit$surv, time = sfit$time, 
#'   majority = c(rep(0, sfit$strata[1]), rep(1, sfit$strata[2]))
#' )
#' 
#' rf_dt <- rbind(
#'   data.table(surv = colMeans(fsurv$srv$srvx[dat$majority == 0, ]), 
#'              time = fsurv$time_interest, majority = 0),
#'   data.table(surv = colMeans(fsurv$srv$srvx[dat$majority == 1, ]), 
#'              time = fsurv$time_interest, majority = 1)
#' )
#' 
#' rfw_dt <- rbind(
#'   data.table(surv = colMeans(fsurv$srv$srvx[dat_run$majority == 0, ]), 
#'              time = fsurv2$time_interest, majority = 0),
#'   data.table(surv = colMeans(fsurv$srv$srvx[dat_run$majority == 1, ]), 
#'              time = fsurv2$time_interest, majority = 1)
#' )
#' 
#' ggplot(rbind(km_dt[, method := "KM"], rf_dt[, method := "RF"],
#'              rfw_dt[, method := "RF-rw"]), 
#'        aes(x=time, y=surv, color = factor(majority), linetype = factor(method))) + 
#'   geom_line(linewidth=1) + theme_bw() +
#'   coord_cartesian(xlim=c(0, 500), ylim = c(0.89, 0.95)) +
#'   scale_color_discrete(name = "Majority") +
#'   scale_linetype_manual(name = "Inference Method", 
#'                         values = c("solid", "dashed", "dotted")) +
#'   theme(legend.position = "inside", legend.position.inside = c(0.65, 0.65),
#'         legend.box.background = element_rect())
#' 
#' 
#' ggplot(
#'   rbind(km_dt[, method := "KM"], rf_dt[, method := "RF"], 
#'         rfw_dt[, method := "RF-rw"])[, diff(surv), by = c("time", "method")],
#'   aes(x = time, y = V1, color = factor(method))
#' ) + geom_line(linewidth=1) +
#'   scale_y_continuous(labels=scales::percent) + theme_bw() +
#'   scale_color_discrete(name="Method") + theme(legend.position="bottom")
#' 
#' # fit the forest to minority data only
#' x0marg <- rfsrc(Surv(event_time, event) ~ age  + sex + irsad + anz_cmb + 
#'                   frailty + apache_iii_diag + elective + apache_iii_rod,
#'                 data = dat[majority == 0])
#' 
#' plot(survfit(Surv(event_time, event) ~ 1, data = dat[majority == 0]),
#'      ylim=c(0.9, 1))
#' lines(x0marg$time.interest, colMeans(x0marg$survival.oob), col="blue")
#' lines(fsurv2$time_interest, colMeans(fsurv2$srv$srvx[dat_run$majority==0, ]), 
#'       col="red")