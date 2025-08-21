
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

data(follic, package = "randomForestSRC")
follic$majority <- rbinom(nrow(follic), 1, prob = 0.5)
follic$rt <- NULL
follic$ch <- as.integer(follic$ch == "Y")
follic <- as.data.table(follic)
follic <- follic[1:100]

# standard fairness model
X <- "majority"
Z <- "age" 
W <- c("clinstg", "hgb", "ch")
time_var <- "time"
event_var <- "status"

#' * classical approach * 
dat_cl <- follic
dat_cl$status <- as.integer(dat_cl$status > 0)
fsrv_cl <- fair_surv(
  dat_cl, X = X, Z = Z, W = W, time_var = time_var, event_var = event_var, 
  nboot = 3
)

str(fsrv_cl)

autoplot(fsrv_cl, scale = "surv")
autoplot(fsrv_cl, scale = "mst")
autoplot(fsrv_cl, scale = "chf")
autoplot(fsrv_cl, scale = "chf-ratio")

#' * competing risks * 
dat_cr <- follic
fsrv_cr <- fair_surv(
  dat_cr, X = X, Z = Z, W = W, time_var = time_var, event_var = event_var, 
  nboot = 3
)

autoplot(fsrv_cr, scale = "surv")
autoplot(fsrv_cr, scale = "cif")

#' * copula model * 
dat_cp <- follic
fsrv_cp <- fair_surv(
  dat_cp, X = X, Z = Z, W = W, time_var = time_var, event_var = event_var, 
  nboot = 3, copula = "clayton", tau = c(0.3, 0.8)
)

autoplot(fsrv_cp, scale = "surv")

