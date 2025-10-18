
# get Kaplan-Meier out-of-bag (5-fold cross-validation)
km_oob <- function(data) {
  
  K <- 5
  surv_oob <- rep(NA, nrow(data))
  folds <- rep_len(seq_len(K), nrow(data))
  for (fold in seq_len(K)) {
    
    fold_idx <- folds == fold
    fit <- survfit(Surv(event_time, event) ~ 1, data = data[!fold_idx])
    
    surv_oob[fold_idx] <- summary(fit, times = data[fold_idx]$event_time, 
                                  extend = TRUE)$surv
  }
  surv_oob
  
  iso <- isoreg(data$event_time, -surv_oob)
  unique(data.table(time = iso$x[iso$ord], surv = -iso$yf, method="KM OOB"))
}
