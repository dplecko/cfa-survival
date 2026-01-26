
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

out <- "death"
dat <- load_data("aics", outcome = out)
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")

datx1 <- dat[majority == 1]

dat_run <- datx1[sample.int(nrow(datx1), size = 10^4)]

dat_run

rhs <- c(X, Z, W)
frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                          paste(rhs, collapse = "+")))

dt_rft <- c()
for (nodesize in c(15, 50, 100, 500, 10^4)) {
  
  rfs <- rfsrc(frml, data = dat_run, nodesize = nodesize)
  dt_rf <- data.table(time = rfs$time.interest, surv = colMeans(rfs$survival.oob),
                      method = "RF", nodesize = nodesize)
  dt_rft <- rbind(dt_rft, dt_rf)
  
}

sfit <- survfit(Surv(event_time, event) ~ 1, data = dat_run)
ridx <- sfit$n.risk > 50
dt_km <- data.table(surv = sfit$surv[ridx], time = sfit$time[ridx], 
                    method = "KM", nodesize = 0)

ggplot(rbind(dt_rft, dt_km), aes(x = time, y=surv, color = factor(nodesize))) +
  geom_line() + theme_bw()
