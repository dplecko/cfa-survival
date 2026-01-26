
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

chf_01 <- function(data, X, time_var, event_var, rhs, method, time_interest,
                   balance_groups, split_forest, ...) {
  
  chf_rfs_cf(data, X, time_var, event_var, rhs, time_interest,
             balance_groups, split_forest, ...)
}

fair_surv <- function(data, X, Z, W, time_var, event_var, rhs=".", 
                      x0 = 0, x1 = 1, 
                      method = c("rfs-cf", "ranger", "ranger-cf", "rfs"),
                      time_interest = NULL,
                      norm_method = "adapt", balance_groups = FALSE,
                      split_forest = FALSE,
                      nboot = 1, copula = NULL, tau_grid = NULL, ...) {
  
  method <- match.arg(method, c("rfs-cf", "ranger", "ranger-cf", "rfs"))
  
  if (nboot > 1L) {
    
    fsrvi <- list()
    pb <- txtProgressBar(min = 0, max = nboot, style = 3)
    for (i in seq_len(nboot)) {
      
      boot_idx <- if (i == 1) TRUE else sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_idx, ]
      
      fsrvi[[i]] <- fair_surv(data = boot_data, X = X, Z = Z, W = W, 
                              time_var = time_var, event_var = event_var, 
                              rhs = rhs, x0 = x0, x1 = x1, 
                              method = method, time_interest = time_interest, 
                              norm_method = norm_method, 
                              balance_groups = balance_groups,
                              split_forest = split_forest,
                              nboot = if (i == 1) 1 else -1, 
                              copula = copula, tau_grid = tau_grid, ...)
      
      if (i == 1) {
        
        time_interest <- fsrvi[[i]]$time_interest
        fsrv <- fsrvi[[i]]
      }
      
      fsrvi[[i]] <- fsrvi[[i]]$measures
      fsrvi[[i]][, boot := i]
      
      setTxtProgressBar(pb, i)
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
  
  # determine rhs
  reg_vars <- c(Z, W)
  if (!split_forest) reg_vars <- c(X, reg_vars)
  rhs <- paste(reg_vars, collapse = "+")
  
  # fit the survival object and obtain counterfactual (OOB) CHF functions
  c(chf, chfx0, chfx1, srv, srvx0, srvx1, cif, cifx0, cifx1, time_interest) %<-% 
    chf_01(data[, c(X, Z, W, time_var, event_var), with=FALSE],
           X, time_var, event_var, rhs, method, time_interest, balance_groups,
           split_forest, ...)
  
  true_nuiss <- FALSE
  if (true_nuiss) {
    
    x1_ind <- data[[X]] == 1
    s_true <- S_T_potential_curves_from_gen(g, tgrid)
    
    srvx0 <- s_true[["S_x0_wx0"]]
    srvx0[x1_ind, ] <- s_true[["S_x0_wx1"]][x1_ind, ]
    
    srvx1 <- s_true[["S_x1_wx0"]]
    srvx1[x1_ind, ] <- s_true[["S_x1_wx1"]][x1_ind, ]
    
    srv <- srvx0
    srv[x1_ind, ] <- srvx1[x1_ind, ]
  }
  
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
              "Considering running with `balance_groups = TRUE` or ",
              "`split_forest = TRUE`.")
    }
  }
  
  if (mean(data[[X]]) > 0.9 & !balance_groups & nboot != -1)
    tv_test(srv, data[, c(X, time_var, event_var), with=FALSE], X, 
            time_var, event_var, time_interest)
  
  # get the propensity scores - regress X on Z, Z+W
  px_zw <- cv_xgb(data[, c(Z, W), with=FALSE], data[[X]])
  px_z  <- cv_xgb(data[, c(Z),    with=FALSE], data[[X]])
  px    <- mean(data[[X]])
  
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
        matx0 <- cif_copula(cifx0, copula, tau)
        matx1 <- cif_copula(cifx1, copula, tau)
      } else {
        
        mat <- scale_transform(scale, chf, srv, cif, time_interest)
        matx0 <- scale_transform(scale, chfx0, srvx0, cifx0, time_interest)
        matx1 <- scale_transform(scale, chfx1, srvx1, cifx1, time_interest)
      }
      
      if (is.element(scale, c("chf", "chf-ratio", "cif"))) {
        
        nslices <- nevents
      } else nslices <- 1
      
      for (eventi in seq_len(nslices)) {
        
        for (i in seq_along(time_interest)) {
          
          # extract fitted values
          if (nslices > 1) {
            
            pout <- mat[, i, eventi]
            poutx0 <- matx0[, i, eventi]
            poutx1 <- matx1[, i, eventi]
          } else {

            pout <- mat[, i]
            poutx0 <- matx0[, i]
            poutx1 <- matx1[, i]
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
  
  res <- melt(res_diff, id.vars = by_cols, variable.name = "effect",
              variable.factor = FALSE)
  if (!is.null(res_ratio)) res <-
    rbind(res, melt(res_ratio, id.vars = by_cols, variable.name = "effect",
                    variable.factor = FALSE))
  
  structure(
    list(
      measures = res, data = data,
      is_cr = is_cr, is_sens = is_sens, copula = copula, tau_grid = tau_grid,
      x0 = x0, x1 = x1, X = X, W = W, Z = Z, time_var = time_var, 
      event_var = event_var, cl = match.call(),
      method = method, time_interest = time_interest, nboot = nboot,
      balance_groups = balance_groups,
      chf = list(chfx = chf, chfx0 = chfx0, chfx1 = chfx1),
      srv = list(srvx = srv, srvx0 = srvx0, srvx1 = srvx1),
      cif = list(cifx = cif, cifx0 = cifx0, cifx1 = cifx1),
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
    inv_gen <- function(t) -1/theta * log ( 1 + exp(-t) * (exp(-theta) - 1) )
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
    
    # H(t, c) is the joint survival
    
    # Step 1: update chat first (this gives a lower bound on chat)
    H_t_tp <- srv[, i] + d_cif2[, i-1]
    chat_ilwr <- inv_gen(gen(H_t_tp) - gen(shat[, i-1]))
    shat_iupr <- inv_gen(
      gen(srv[, i]) - gen(srv[, i-1]) + gen(chat[, i-1]) - gen(chat_ilwr) +
        gen(shat[, i-1])
    )
    
    # Step 2: update shat first (this gives a lower bound on shat)
    H_tp_t <- srv[, i] + d_cif1[, i-1]
    shat_ilwr <- inv_gen(gen(H_tp_t) - gen(chat[, i-1]))
    chat_iupr <- inv_gen(
      gen(srv[, i]) - gen(srv[, i-1]) + gen(shat[, i-1]) - gen(shat_ilwr) +
        gen(chat[, i-1])
    )
    
      
    chat[, i] <- 1/2 * (chat_ilwr + chat_iupr)
    shat[, i] <- 1/2 * (shat_ilwr + shat_iupr)
  }
  
  shat[, -1]
}

cif_copula_single <- function(cif1, cif2, copula, tau) {
  # cif1, cif2: numeric vectors over the same time grid (no t=0 included)
  stopifnot(is.numeric(cif1), is.numeric(cif2), length(cif1) == length(cif2))
  
  theta <- tau_to_theta(tau, copula)
  gen <- copula_gen(theta, copula)
  inv_gen <- copula_gen(theta, copula, inv = TRUE)
  
  cif1 <- c(0, cif1)
  cif2 <- c(0, cif2)
  
  srv <- 1 - cif1 - cif2
  d1 <- diff(cif1)
  d2 <- diff(cif2)
  
  m <- length(srv)
  shat <- chat <- rep(1, m)
  
  for (i in 2:m) {
    
    # Step 1: update chat first
    H_t_tp <- srv[i] + d2[i-1]
    chat_ilwr <- inv_gen(gen(H_t_tp) - gen(shat[i-1]))
    shat_iupr <- inv_gen(
      gen(srv[i]) - gen(srv[i-1]) + gen(chat[i-1]) - gen(chat_ilwr) + gen(shat[i-1])
    )
    
    # Step 2: update shat first
    H_tp_t <- srv[i] + d1[i-1]
    shat_ilwr <- inv_gen(gen(H_tp_t) - gen(chat[i-1]))
    chat_iupr <- inv_gen(
      gen(srv[i]) - gen(srv[i-1]) + gen(shat[i-1]) - gen(shat_ilwr) + gen(chat[i-1])
    )
    
    if (any(is.nan(c(chat_ilwr, chat_iupr, shat_ilwr, shat_iupr)))) browser()
    
    # Step 3: take the midpoint of the two bounds
    chat[i] <- 0.5 * (chat_ilwr + chat_iupr)
    shat[i] <- 0.5 * (shat_ilwr + shat_iupr)
  }
  
  shat[-1]
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

#' * new fair_surv code -- modular IPW/model-based/AIPW *
fs_fit <- function(data, X, Z, W, time_var, event_var,
                   time_interest = NULL,
                   balance_groups = FALSE, split_forest = FALSE, ...) {
  # propensities (exactly as before)
  e_zw <- cv_xgb(data[, c(Z, W), with = FALSE], data[[X]])
  e_z  <- cv_xgb(data[, c(Z),    with = FALSE], data[[X]])
  p1   <- mean(data[[X]])
  
  # outcome S/CIF via your rfs learner
  rhs_out <- paste(c(if (!split_forest) X, Z, W), collapse = "+")
  out <- chf_rfs_cf(
    data = data[, c(X, Z, W, time_var, event_var), with = FALSE],
    X = X, time_var = time_var, event_var = event_var,
    rhs = rhs_out, time_interest = time_interest,
    balance_groups = balance_groups, split_forest = split_forest, ...
  )
  
  # censoring G(t): flip event to "censoring event" (1 if censored)
  dG <- copy(data)
  dG[[event_var]] <- as.integer(dG[[event_var]] == 0L)
  rhs_cen <- paste(c(if (!split_forest) X, Z, W), collapse = "+")
  cen <- chf_rfs_cf(
    data = dG[, c(X, Z, W, time_var, event_var), with = FALSE],
    X = X, time_var = time_var, event_var = event_var,
    rhs = rhs_cen, time_interest = out$time_interest,
    balance_groups = balance_groups, split_forest = split_forest, ...
  )
  
  list(
    data = data,
    X = X, Z = Z, W = W, time_var = time_var, event_var = event_var,
    grid = list(t = out$time_interest, J = if (is.null(out$cif)) 1L else dim(out$cif)[3]),
    ps = list(e_z = e_z, e_zw = e_zw, p1 = p1, meta = list(engine = "xgb")),
    S = list(
      S_x0 = out$srvx0, S_x1 = out$srvx1,
      CIF_x0 = out$cifx0, CIF_x1 = out$cifx1,
      meta = list(engine = "rfsrc")
    ),
    G = list(
      G_x0 = cen$srvx0, G_x1 = cen$srvx1,
      meta = list(engine = "rfsrc")
    )
  )
}

fs_estimate <- function(pack, method = c("model","ipw","aipw")) {
  method <- match.arg(method)
  switch(method,
         model = fs_estimate_model(pack),
         ipw   = fs_estimate_ipw(pack),
         aipw  = fs_estimate_aipw(pack)
  )
}

fair_surv_v2 <- function(data, X, Z, W, time_var, event_var,
                         time_interest = NULL,
                         balance_groups = FALSE, split_forest = FALSE,
                         method = c("model","ipw","aipw"), ...) {
  
  pack <- fs_fit(data, X, Z, W, time_var, event_var,
                 time_interest, balance_groups, split_forest, ...)
  est  <- fs_estimate(pack, method = match.arg(method))
  est
}
