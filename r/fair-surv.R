
match_grids <- function(a, b) {
  
  res <- integer(length(a))
  bi <- 1
  for (ai in seq_along(a)) {
    while (bi <= length(b) && b[bi] < a[ai]) {
      bi <- bi + 1
    }
    if (bi > length(b)) {
      res[ai] <- length(b)  # a[i] > max(b)
    } else if (b[bi] < a[ai]) {
      res[ai] <- NA_integer_
    } else if (bi == 1 && b[bi] > a[ai]) {
      res[ai] <- NA_integer_  # a[i] < min(b)
    } else {
      res[ai] <- bi
    }
  }
  
  res
}

match_grids_lwr <- function(a, b) {
  
  res <- rep(NA_integer_, length(a))
  acurr <- bcurr <- 1
  while (acurr <= length(a) & bcurr <= length(b)) {
    
    # assign index if possible
    if (a[acurr] >= b[bcurr]) {
      
      res[acurr] <- bcurr
      acurr <- acurr + 1
    }
    
    # move b to the left as much as possible
    while(bcurr+1 <= length(b) && acurr <= length(a) && b[bcurr+1] <= a[acurr]) {
      bcurr <- bcurr + 1
    }
  }
  
  res
}

cv_xgb <- function(df, y, weights = NULL, ...) {
  
  if (is.character(as.matrix(df))) browser()
  dtrain <- xgb.DMatrix(data = as.matrix(df), label = y, weight = weights)
  
  binary <- all(y %in% c(0, 1))
  if (binary) {
    
    params <- list(objective = "binary:logistic", eval_metric = "logloss")
  } else {
    
    params <- list(objective = "reg:squarederror", eval_metric = "rmse")
  }
  
  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    prediction = TRUE,
    verbose = FALSE, ...
  )
  
  return(cv$pred)
}

chf_ranger <- function(data, X, time_var, event_var, rhs, time_interest,
                       balance_groups) {
  
  data_x0 <- copy(data)
  data_x1 <- copy(data)
  data_x0[[X]] <- 0
  data_x1[[X]] <- 1
  
  frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  obj <- ranger(frml, data, keep.inbag = TRUE, time.interest = time_interest)
  
  pred <- predict(obj, data = data, predict.all = TRUE)$chf
  pred_x0 <- predict(obj, data = data_x0, predict.all = TRUE)$chf
  pred_x1 <- predict(obj, data = data_x1, predict.all = TRUE)$chf
  pred <- aperm(pred, perm = c(1, 3, 2))
  pred_x0 <- aperm(pred_x0, perm = c(1, 3, 2))
  pred_x1 <- aperm(pred_x1, perm = c(1, 3, 2))
  
  oob_ind <- do.call(cbind, lapply(obj$inbag.counts, function(x) x == 0))
  oob3 <- array(oob_ind, dim = dim(pred_x0))
  
  chf <- colSums(aperm(pred * oob3, c(2, 1, 3))) / rowSums(oob_ind)
  chf_x0 <- colSums(aperm(pred_x0 * oob3, c(2, 1, 3))) / rowSums(oob_ind)
  chf_x1 <- colSums(aperm(pred_x1 * oob3, c(2, 1, 3))) / rowSums(oob_ind)
  
  list(
    chf = chf,
    chf_x0 = chf_x0,
    chf_x1 = chf_x1,
    time_interest = obj$unique.death.times
  )
}

chf_ranger_cf  <- function(data, X, time_var, event_var, rhs, time_interest, 
                           balance_groups, K = 5) {

  n <- nrow(data)
  idx <- sample(rep(1:K, length.out = n))
  
  chf <- chf_x0 <- chf_x1 <- vector("list", K)
  if (is.null(time_interest)) time_interest <- 150 # use 150 time points
  ind <- NULL
  
  for (k in 1:K) {
    trn <- data[idx != k]
    val <- data[idx == k]
    ind <- c(ind, which(idx == k))
    
    val_x0 <- copy(val)
    val_x1 <- copy(val)
    val_x0[[X]] <- 0
    val_x1[[X]] <- 1
    
    frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
    obj <- ranger(frml, data = trn, time.interest = time_interest)
    
    p <- predict(obj, data = val)
    p0 <- predict(obj, data = val_x0)
    p1 <- predict(obj, data = val_x1)
    
    chf[[k]] <- p$chf
    chf_x0[[k]] <- p0$chf
    chf_x1[[k]] <- p1$chf
    
    if (length(time_interest) == 1) time_interest <- obj$unique.death.times
  }
  
  list(
    chf = do.call(rbind, chf)[order(ind), ],
    chf_x0 = do.call(rbind, chf_x0)[order(ind), ],
    chf_x1 = do.call(rbind, chf_x1)[order(ind), ],
    time_interest = time_interest
  )
}

chf_rfs <- function(data, X, time_var, event_var, rhs, time_interest,
                    balance_groups) {
  
  frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  obj <- rfsrc(frml, data, ntree = 500, ntime = time_interest,
               nodesize = 5, nsplit = 50, membership = TRUE)
  
  data_x0 <- data_x1 <- copy(data)
  data_x0[[X]] <- 0
  data_x1[[X]] <- 1
  chf_x0 <- chf_x1 <- 0
  
  for (i in seq_len(obj$ntree)) {
    chf_x0 <- chf_x0 + predict(obj, newdata = data_x0, get.tree = i)$chf * 
      (1 - obj$inbag[, i])
    chf_x1 <- chf_x1 + predict(obj, newdata = data_x1, get.tree = i)$chf * 
      (1 - obj$inbag[, i])
  }
  
  chf_x0 <- chf_x0 / rowSums(1 - obj$inbag)
  chf_x1 <- chf_x1 / rowSums(1 - obj$inbag)
  
  list(chf = obj$chf.oob, chf_x0 = chf_x0, chf_x1 = chf_x1, 
       time_interest = obj$time.interest)
}

chf_rfs_cf <- function(data, X, time_var, event_var, rhs, time_interest,
                       balance_groups, K = 5) {
  n <- nrow(data)
  idx <- sample(rep(1:K, length.out = n))
  
  
  if (balance_groups) {
    
    wt0 <- 1 / mean(data[[X]] == 0)
    wt1 <- 1 / mean(data[[X]] == 1) 
  }
  
  chf <- chf_x0 <- chf_x1 <- srv <- srv_x0 <- srv_x1 <- 
    cif <- cif_x0 <- cif_x1 <- vector("list", K)
  if (is.null(time_interest)) time_interest <- 150 # use 150 time points
  ind <- NULL
  
  # check if any competing risk events are considered
  is_cr <- if (!all(data[[event_var]] %in% c(0, 1))) TRUE else FALSE
  
  for (k in 1:K) {
    trn <- data[idx != k]
    val <- data[idx == k]
    ind <- c(ind, which(idx == k))
    
    val_x0 <- copy(val)
    val_x1 <- copy(val)
    val_x0[[X]] <- 0
    val_x1[[X]] <- 1
    
    frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
    
    # training weights
    gwt <- if (balance_groups) ifelse(trn[[X]], wt1, wt0) else rep(1, nrow(trn))
    gsize <- if (balance_groups) 2 * min(table(trn[[X]])) else nrow(trn)
    obj <- rfsrc(frml, data = trn, ntime = time_interest, case.wt = gwt,
                 samptype = "swr", sampsize = gsize)
    
    preds <- predict(obj, newdata = val)
    preds_x0 <- predict(obj, newdata = val_x0)
    preds_x1 <- predict(obj, newdata = val_x1)
    
    chf[[k]] <- preds$chf
    chf_x0[[k]] <- preds_x0$chf
    chf_x1[[k]] <- preds_x1$chf
    
    if (is_cr) {
      
      cif[[k]] <- preds$cif
      cif_x0[[k]] <- preds_x0$cif
      cif_x1[[k]] <- preds_x1$cif
      
      # for competing risks, overall survival S(t) = 1 - \sum CIF_j(t)
      srv[[k]] <- 1 - apply(cif[[k]], c(1, 2), sum)
      srv_x0[[k]] <- 1 - apply(cif_x0[[k]], c(1, 2), sum)
      srv_x1[[k]] <- 1 - apply(cif_x1[[k]], c(1, 2), sum)
    } else {
      
      srv[[k]] <- preds$survival
      srv_x0[[k]] <- preds_x0$survival
      srv_x1[[k]] <- preds_x1$survival
    }
    
    if (length(time_interest) == 1) {
      
      time_interest <- obj$time.interest
    } else {
      
      grid_ord <- match_grids(time_interest, c(-Inf, obj$time.interest))
      
      if (is_cr) {
        
        zeros_dim <- dim(chf[[k]])
        zeros_dim[2] <- 1
        zeros <- array(0, dim = zeros_dim)
        ones <- array(1, dim = zeros_dim)
        
        chf[[k]] <- abind(zeros, chf[[k]], along = 2)[, grid_ord, ]
        chf_x0[[k]] <- abind(zeros, chf_x0[[k]], along = 2)[, grid_ord, ]
        chf_x1[[k]] <- abind(zeros, chf_x1[[k]], along = 2)[, grid_ord, ]
        
        srv[[k]] <- cbind(1, srv[[k]])[, grid_ord]
        srv_x0[[k]] <- cbind(1, srv_x0[[k]])[, grid_ord]
        srv_x1[[k]] <- cbind(1, srv_x1[[k]])[, grid_ord]
        
        cif[[k]] <- abind(ones, cif[[k]], along = 2)[, grid_ord, ]
        cif_x0[[k]] <- abind(ones, cif_x0[[k]], along = 2)[, grid_ord, ]
        cif_x1[[k]] <- abind(ones, cif_x1[[k]], along = 2)[, grid_ord, ]
        
      } else {
        
        chf[[k]] <- cbind(0, chf[[k]])[, grid_ord]
        chf_x0[[k]] <- cbind(0, chf_x0[[k]])[, grid_ord]
        chf_x1[[k]] <- cbind(0, chf_x1[[k]])[, grid_ord]
        
        srv[[k]] <- cbind(1, srv[[k]])[, grid_ord]
        srv_x0[[k]] <- cbind(1, srv_x0[[k]])[, grid_ord]
        srv_x1[[k]] <- cbind(1, srv_x1[[k]])[, grid_ord]
      }
    }
  }
  
  if (is_cr) {
    
    ret <- list(
      chf = do.call(abind, args = list(chf, along = 1))[order(ind), ,],
      chf_x0 = do.call(abind, args = list(chf_x0, along = 1))[order(ind), ,],
      chf_x1 = do.call(abind, args = list(chf_x1, along = 1))[order(ind), ,],
      srv = do.call(rbind, srv)[order(ind), ],
      srv_x0 = do.call(rbind, srv_x0)[order(ind), ],
      srv_x1 = do.call(rbind, srv_x1)[order(ind), ],
      cif = do.call(abind, args = list(cif, along = 1))[order(ind), ,],
      cif_x0 = do.call(abind, args = list(cif_x0, along = 1))[order(ind), ,],
      cif_x1 = do.call(abind, args = list(cif_x1, along = 1))[order(ind), ,],
      time_interest = time_interest
    )
  } else {
    
    ret <- list(
      chf = do.call(rbind, chf)[order(ind), ],
      chf_x0 = do.call(rbind, chf_x0)[order(ind), ],
      chf_x1 = do.call(rbind, chf_x1)[order(ind), ],
      srv = do.call(rbind, srv)[order(ind), ],
      srv_x0 = do.call(rbind, srv_x0)[order(ind), ],
      srv_x1 = do.call(rbind, srv_x1)[order(ind), ],
      cif = NULL, cif_x0 = NULL, cif_x1 = NULL,
      time_interest = time_interest
    )
  }
  
  ret
}

chf_01 <- function(data, X, time_var, event_var, rhs, method, time_interest,
                   balance_groups) {
  
  if (method == "ranger") 
    return(chf_ranger(data, X, time_var, event_var, rhs, time_interest,
                      balance_groups))
  if (method == "ranger-cf") 
    return(chf_ranger_cf(data, X, time_var, event_var, rhs, time_interest,
                         balance_groups))
  if (method == "rfs") 
    return(chf_rfs(data, X, time_var, event_var, rhs, time_interest,
                   balance_groups))
  if (method == "rfs-cf") 
    return(chf_rfs_cf(data, X, time_var, event_var, rhs, time_interest,
                      balance_groups))
}

fair_surv <- function(data, X, Z, W, time_var, event_var, rhs=".", 
                      x0 = 0, x1 = 1, 
                      method = c("rfs-cf", "ranger", "ranger-cf", "rfs"),
                      time_interest = NULL,
                      norm_method = "adapt", balance_groups = FALSE, 
                      nboot = 1, copula = NULL, tau_grid = NULL, ...) {
  
  method <- match.arg(method, c("rfs-cf", "ranger", "ranger-cf", "rfs"))
  
  if (nboot > 1L) {
    
    fsrvi <- list()
    for (i in seq_len(nboot)) {
      
      boot_idx <- if (i == 1) TRUE else sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_idx, ]
      
      fsrvi[[i]] <- fair_surv(data = boot_data, X = X, Z = Z, W = W, 
                              time_var = time_var, event_var = event_var, 
                              rhs = rhs, x0 = x0, x1 = x1, 
                              method = method, time_interest = time_interest, 
                              norm_method = norm_method, 
                              balance_groups = balance_groups,
                              nboot = if (i == 1) 1 else -1, 
                              copula = copula, tau_grid = tau_grid, ...)
      
      if (i == 1) {
        
        time_interest <- fsrvi[[i]]$time_interest
        fsrv <- fsrvi[[i]]
      }
      
      fsrvi[[i]] <- fsrvi[[i]]$measures
      fsrvi[[i]][, boot := i]
    }
    
    fsrv_meas <- do.call(rbind, fsrvi)
    by_cols <- c("time_interest", "effect", "scale", "event")
    if (!is.null(tau_grid) & !is.null(copula)) 
      by_cols <- c(by_cols, "tau", "copula")
    
    fsrv_meas <- fsrv_meas[, list(value = value[boot == 1], sd = sd(value)), 
                           by = by_cols]
    fsrv$measures <- fsrv_meas
    return(fsrv)
  }
  
  # determine the number of competing events
  nevents <- max(data[[event_var]])
  is_sens <- !is.null(copula) || !is.null(tau_grid)
  is_cr <- nevents > 1 & !is_sens
  
  # fit the survival object and obtain counterfactual (OOB) CHF functions
  c(chf, chf_x0, chf_x1, srv, srv_x0, srv_x1, cif, cif_x0, cif_x1, time_interest) %<-% 
    chf_01(data[, c(X, Z, W, time_var, event_var), with=FALSE],
           X, time_var, event_var, rhs, method, time_interest, balance_groups)
  
  tv_test <- function(srv, data, X, time_var, event_var, time_interest) {
    
    idx1 <- data[[X]] == 1
    margx1 <- colMeans(srv[idx1, ])
    margx0 <- colMeans(srv[!idx1, ])
    tv_rf <- margx1 - margx0
    
    frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ 1"))
    fit1 <- survfit(frml, data = data[get(X) == 1])
    sd1 <- (fit1$surv - fit1$lower) / 1.96
    fit0 <- survfit(frml, data = data[get(X) == 0])
    sd0 <- (fit0$surv - fit0$lower) / 1.96
    
    match0 <- match_grids_lwr(time_interest, fit1$time)
    match1 <- match_grids_lwr(time_interest, fit0$time)
    
    tv_km <- (fit1$surv[match1] - fit0$surv[match0])
    sd_tv <- sqrt(sd1[match1]^2 + sd0[match0]^2)
    
    
    viol <- any(tv_rf < tv_km - 1.96 * sd_tv | tv_rf > tv_km + 1.96 * sd_tv)
    if (viol) {
      
      message("TV estimate with RF different from Kaplan-Meier.\n",
              "This may be due to imbalance of groups according to X.\n",
              "Considering running with `balance_groups = TRUE`.")
    }
  }
  
  if (mean(data[[X]]) > 0.9 & !balance_groups & nboot != -1)
    tv_test(srv, data[, c(X, time_var, event_var), with=FALSE], X, 
            time_var, event_var, time_interest)
  
  # get the propensity scores - regress X on Z, Z+W
  px_zw <- cv_xgb(data[, c(Z, W), with=FALSE], data[[X]])
  px_z <- cv_xgb(data[, c(Z), with=FALSE], data[[X]])
  px <- mean(data[[X]])
  
  wgh_sum <- function(x, wgh) sum(x * wgh) / sum(wgh)
  
  compute_po <- function(fx = 0, wx = 0, zx = 0, thresh = FALSE, 
                         norm_method = "adapt") {
    
    if (is.na(zx)) zx <- -1
    
    if (wx == 0) wgh <- (1 - px_zw) / (1 - px_z) else wgh <- px_zw / px_z
    if (zx == 0) {
      wgh <- wgh * (1 - px_zw) / (1 - px)
    } else if (zx == 1) wgh <- wgh * px_z / px
    
    if (fx == 0) po_samp <- poutx0 else po_samp <- poutx1
    if (norm_method == "adapt") norm_const <- sum(wgh) else 
      norm_const <- length(po_samp)
    if (thresh) po_samp <- as.integer(po_samp > thr)
    list(po_samp = po_samp, wgh = wgh, norm_const = norm_const)
  }
  
  eval_po <- function(po) sum(po$po_samp * po$wgh) / po$norm_const
  diff_po <- function(po1, po2) eval_po(po1) - eval_po(po2)
  ratio_po <- function(po1, po2) eval_po(po1) / eval_po(po2)
  
  scale_transform <- function(scale, chf, srv, cif, time_interest) {
    
    if (scale == "surv") {
      
      mat <- srv
    } else if (scale == "mst") {
      
      cumsum(t(t(srv) * diff(c(0, time_interest))))
    } else if (is.element(scale, c("chf", "chf-ratio"))) {
      
      mat <- chf
    } else if (scale == "cif") {
      
      mat <- cif
    }

    return(mat)
  }
  
  res_diff <- res_ratio <- c()
  scales <- c("surv", "mst", "chf", "chf-ratio")
  if (is_cr) scales <- c(scales, "cif")
  if (is_sens) scales <- "surv"
  
  #' * tau loop for the sensitivity case *
  if (is.null(tau_grid)) tau_grid <- 0
  for (tau in tau_grid) {
    
    for (scale in scales) {
      
      if (is_sens) {
        
        # convert CIFs to S_A
        mat <- cif_copula(cif, copula, tau)
        mat_x0 <- cif_copula(cif_x0, copula, tau)
        mat_x1 <- cif_copula(cif_x1, copula, tau)
      } else {
        
        mat <- scale_transform(scale, chf, srv, cif, time_interest)
        mat_x0 <- scale_transform(scale, chf_x0, srv_x0, cif_x0, time_interest)
        mat_x1 <- scale_transform(scale, chf_x1, srv_x1, cif_x1, time_interest)
      }
      
      if (is.element(scale, c("chf", "chf-ratio", "cif"))) {
        
        nslices <- nevents
      } else nslices <- 1
      
      for (eventi in seq_len(nslices)) {
        
        for (i in seq_along(time_interest)) {
          
          # extract fitted values
          if (nslices > 1) {
            
            pout <- mat[, i, eventi]
            poutx0 <- mat_x0[, i, eventi]
            poutx1 <- mat_x1[, i, eventi]
          } else {
            
            pout <- mat[, i]
            poutx0 <- mat_x0[, i]
            poutx1 <- mat_x1[, i]
          }
          
          # get f_{x_1, W_{x_0}}
          fx1wx0 <- compute_po(fx = 1, wx = 0, zx = NA, FALSE, norm_method)
          fx1wx0_x0 <- compute_po(fx = 1, wx = 0, zx = 0, FALSE, norm_method)
          
          # get f_{x_0, W_{x_0}}
          fx0wx0 <- compute_po(fx = 0, wx = 0, zx = NA, FALSE, norm_method)
          fx0wx0_x0 <- compute_po(fx = 0, wx = 0, zx = 0, FALSE, norm_method)
          
          # get f_{x_1, W_{x_1}}
          fx1wx1 <- compute_po(fx = 1, wx = 1, zx = NA, FALSE, norm_method)
          fx1wx1_x0 <- compute_po(fx = 1, wx = 1, zx = 0, FALSE, norm_method)
          
          # get f | x0
          f_x0 <- mean(pout[data[[X]] == 0])
          fx0wx0_x0 <- compute_po(fx = 0, wx = 0, zx = 0, FALSE, norm_method)
          
          # get f | x1
          f_x1 <- mean(pout[data[[X]] == 1])
          fx1wx1_x1 <- compute_po(fx = 1, wx = 1, zx = 1, FALSE, norm_method)
          
          if (scale != "chf-ratio") {
            
            # natural effects
            tv <- f_x1 - f_x0
            nde <- diff_po(fx1wx0, fx0wx0)
            nie <- diff_po(fx1wx0, fx1wx1)
            nse <- diff_po(fx1wx1_x1, fx1wx1) - diff_po(fx0wx0_x0, fx0wx0)
            
            # counterfactual effects
            ctfde <- diff_po(fx1wx0_x0, fx0wx0_x0)
            ctfie <- diff_po(fx1wx0_x0, fx1wx1_x0)
            ctfse <- diff_po(fx1wx1_x0, fx1wx1_x1)
            
            res_add <- data.table(
              tv = tv, nde = nde, nie = nie, nse = nse,
              ctfde = ctfde, ctfie = ctfie, ctfse = ctfse, 
              time_interest = time_interest[i], scale = scale,
              event = eventi
            )
            if (is_sens) {
              
              res_add$tau <- tau
              res_add$copula <- copula
            }
            res_diff <- rbind(res_diff, res_add)
          } else {
            
            ctfdr <- ratio_po(fx1wx0_x0, fx0wx0_x0)
            ctfir <- ratio_po(fx1wx0_x0, fx1wx1_x0)
            ctfsr <- ratio_po(fx1wx1_x0, fx1wx1_x1)
            res_ratio <- rbind(
              res_ratio,
              data.table(
                ctfdr = ctfdr, ctfir = ctfir, ctfsr = ctfsr, 
                time_interest = time_interest[i], scale = scale,
                event = eventi
              )
            )
          }
        }
      }
    }
  }
  
  by_cols <- c("time_interest", "scale", "event")
  if (is_sens) by_cols <- c(by_cols, "tau", "copula")
  
  res <- melt(res_diff, id.vars = by_cols, variable.name = "effect")
  if (!is.null(res_ratio)) res <-
    rbind(res, melt(res_ratio, id.vars = by_cols, variable.name = "effect"))
  
  structure(
    list(
      measures = res, data = data,
      is_cr = is_cr, is_sens = is_sens, copula = copula, tau_grid = tau_grid,
      x0 = x0, x1 = x1, X = X, W = W, Z = Z, time_var = time_var, 
      event_var = event_var, cl = match.call(),
      method = method, time_interest = time_interest, nboot = nboot,
      balance_groups = balance_groups,
      chf = list(chfx = chf, chf_x0 = chf_x0, chf_x1 = chf_x1),
      srv = list(srvx = srv, srv_x0 = srv_x0, srv_x1 = srv_x1),
      cif = list(cifx = cif, cif_x0 = cif_x0, cif_x1 = cif_x1),
      pw = list(px_zw = px_zw, px_z = px_z, px = px)
    ),
    class = "fairsurv"
  )
}


tau_to_theta <- function(tau, copula) {
  
  cid <- switch(copula, clayton = 3, gumbel = 4, frank = 5)
  BiCopTau2Par(cid, tau)
}

copula_gen <- function(theta, copula, inv = FALSE) {
  
  if (copula == "clayton") {
    
    gen <- function(t) 1/theta * (t^(-theta) - 1)
    inv_gen <- function(t) (1 + theta * t)^(-1/theta)
  } else if (copula == "gumbel") {
    
    gen <- function(t) (-log(t))^theta
    inv_gen <- function(t) exp(-t^(1/theta))
  } else if (copula == "frank") {
    
    gen <- function(t) -log ( (exp(-theta * t ) - 1) / (exp(-theta) - 1) )
    inv_gen <- function(t) 1/theta * log ( 1 + exp(-t) * (exp(-theta) - 1) )
  }
  
  if (inv) return(inv_gen) else return(gen)
}

cif_copula <- function(cif, copula, tau) {
  
  theta <- tau_to_theta(tau, copula)
  gen <- copula_gen(theta, copula)
  inv_gen <- copula_gen(theta, copula, inv = TRUE)
  cif1 <- cbind(0, cif[, , 1])
  cif2 <- cbind(0, cif[, , 2])
  srv <- 1 - cif1 - cif2
  d_cif1 <- t(diff(t(cif1)))
  d_cif2 <- t(diff(t(cif2)))
  
  shat <- chat <- array(1, dim = dim(srv))
  
  for (i in seq.int(2, ncol(srv))) {
    
    # # Delta CIF2 = 0 value
    # inv_gen(gen(srv[, i]) - gen(srv[, i-1]) + gen(shat[, i-1])) # shat
    # chat[, i-1] # chat
    # 
    # # Delta CIF1 = 0 value
    # shat[, i-1] # shat
    # inv_gen(gen(srv[, i]) - gen(srv[, i-1]) + gen(chat[, i-1])) # chat value
    
    # both Delta CIF != 0
    sjnt_t_tp <- srv[, i] + d_cif2[, i-1]
    chat[, i] <- inv_gen(gen(sjnt_t_tp) - gen(shat[, i-1])) # chat
    shat[, i] <- inv_gen(
      gen(srv[, i]) - gen(srv[, i-1]) + gen(chat[, i-1]) - gen(chat[, i]) +
        gen(shat[, i-1])
    )
  }
  
  shat[, -1]
}

autoplot.fairsurv <- function(object, 
                              scale = c("surv", "mst", "chf", "chf-ratio", "cif"), 
                              lvl = c("xspecific", "population"), ...) {
  
  scale_val <- match.arg(scale, c("surv", "mst", "chf", "chf-ratio", "cif"))
  lvl <- match.arg(lvl, c("xspecific", "population"))
  ttl <- switch (
    scale_val,
    surv = "Survival", mst = "Mean Survival Time", chf = "Cumulative Hazard",
    `chf-ratio` = "Cumulative Hazard Ratio", cif = "Cumulative Incidence Function"
  )
  
  if (scale_val == "chf-ratio") {
    
    meas <- c("ctfdr", "ctfir", "ctfsr")
  } else if (lvl == "xspecific") {
    
    meas <- c("ctfde", "ctfie", "ctfse")
  } else meas <- c("nde", "nie", "nse")
  
  alpha <- 0.05 # fixed for now
  width <- qnorm(1 - alpha / 2)
  effect_lab <- if (scale_val == "chf-ratio") "Ratio" else "Effect"
  
  plt_dat <- object$measures[scale == scale_val & effect %in% meas]
  
  if (length(unique(plt_dat$event)) > 1) 
    plt_dat$event <- paste("Event =", plt_dat$event)
  
  if (length(unique(plt_dat$tau)) > 1) 
    plt_dat$tau <- paste("Tau =", plt_dat$tau)
  
  p <- ggplot(plt_dat, aes(x = time_interest, y = value, color = effect)) +
    geom_line(linewidth = 1) + theme_bw() +
    xlab("Time") + ylab(effect_lab) +
    scale_color_discrete(
      name = effect_lab, labels = c("Direct", "Indirect", "Spurious")
    ) +
    scale_fill_discrete(
      name = effect_lab, labels = c("Direct", "Indirect", "Spurious")
    ) + ggtitle(ttl)
  
  if (length(unique(plt_dat$event)) > 1) p <- p + facet_wrap(~event)
  if (length(unique(plt_dat$tau)) > 1) p <- p + facet_wrap(~tau)
  
  if (is.element("sd", names(object$measures))) {
    
    p <- p + geom_ribbon(aes(ymin = value - width * sd, 
                             ymax = value + width * sd,
                             fill = effect),
                         alpha = 0.3, linewidth = 0)
  }
  
  if (scale_val == "chf-ratio" & length(plt_dat$event) == 1) {
    cox_data <- object$data[, c(object$X, object$Z, object$W, object$event_var, 
                                object$time_var), 
                            with=FALSE]
    frml <- as.formula(paste0("Surv(", object$time_var, ", ", 
                              object$event_var, ") ~ ."))
    coxmod <- coxph(frml, data = cox_data)
    est <- coxmod$coefficients[[object$X]]
    se <- sqrt(vcov(coxmod)[object$X, object$X])
    ci <- exp(est + c(-1, 1) * width * se)
    p <- p + geom_hline(yintercept = exp(coxmod$coefficients[[object$X]]), 
                        color = "orange", linetype = "dashed") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = ci[1], ymax = ci[2],
               alpha = 0.1, fill = "orange") +
      annotate("text", x = Inf, y = exp(est), label = "Cox-PH Direct CHF Ratio",
               hjust = 1.1, vjust = -0.5, color = "orange")
    
  }
  
  p
}
