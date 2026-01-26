
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# select data source
src <- "aics"
out <- "death" # "dcr", "readm"
balance <- FALSE
splitf <- FALSE

# prepare the data and the SFM
dat <- load_data("aics", outcome = out)
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")

set.seed(2026)

# for testing
dat_run <- dat[sample.int(nrow(dat), size = 1000, replace=FALSE)]

fsurv <- fair_surv(dat_run, X, Z, W, time_var, event_var, 
                   balance_groups = balance, split_forest = splitf)

fsurv_v2m <- fair_surv_v2(dat_run, X, Z, W, time_var, event_var, 
                      balance_groups = balance,
                      split_forest = splitf,
                      method = "model")

fsurv_v2i <- fair_surv_v2(dat_run, X, Z, W, time_var, event_var, 
                      balance_groups = balance,
                      split_forest = splitf,
                      method = "ipw")

fsurv_old <- fsurv$measures[scale == "surv"]
fsurv_old[, method := "model-old"]
setnames(fsurv_old, "time_interest", "time")

ggplot(
  rbind(
    fsurv_old,
    fsurv_v2m$measures[, method := "model-based"]#,
    #fsurv_v2i$measures[, method := "ipw"]
  ), aes(x = time, y = value, color = method)
) + geom_line(linewidth=1.5) + theme_bw() +
  facet_wrap(~ effect)
