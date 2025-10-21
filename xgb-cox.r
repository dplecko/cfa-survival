library(survival); library(xgboost)

set.seed(1)
n <- 400; p <- 5
Z <- matrix(rnorm(n*p), n, p)
b <- c(0.7,-0.5,0.4,0,0); lp0 <- drop(Z %*% b)
Tt <- rexp(n, rate=exp(lp0)); C <- rexp(n, rate=0.3)
y <- pmin(Tt, C); d <- as.integer(Tt <= C)

tr <- sample.int(n, 300); te <- setdiff(seq_len(n), tr)

# xgboost (Cox); encode censoring with negative time
lab_tr <- ifelse(d[tr]==1, y[tr], -y[tr])
dm_tr <- xgb.DMatrix(Z[tr,], label = lab_tr)

par <- list(objective="survival:cox", eval_metric="cox-nloglik",
            eta=0.05, max_depth=3, subsample=0.8, colsample_bytree=0.8)
bst <- xgb.train(par, dm_tr, nrounds=300, verbose=0)

# linear predictors
lp_tr <- predict(bst, Z[tr,])
lp_te <- predict(bst, Z[te,])

# baseline hazard via Cox with offset; Breslow
fit <- coxph(Surv(y[tr], d[tr]) ~ offset(lp_tr), ties="breslow")
bh <- basehaz(fit, centered=FALSE)  # H0(t) step fn

# survival curves S(t|Z) on bh$time grid
S_from_lp <- function(lp, bh) { exp(-outer(bh$hazard, exp(lp))) } # rows=time, cols=obs
S_te <- S_from_lp(lp_te, bh)

# usage examples:
t_idx <- which.max(bh$time >= 1.0)           # S at t≈1
S_at1 <- S_te[t_idx, ]                       # vector length = length(te)
rmst <- function(tt, S) sum(diff(c(0,tt)) * (head(S,-1)+tail(S,-1))/2)
RMST_te <- apply(S_te, 2, rmst, tt = bh$time)  # RMST up to max(bh$time)