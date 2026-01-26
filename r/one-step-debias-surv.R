
cv_xgb_surv <- function(df, y, weights = NULL, ...) {
  
  if (!is.vector(y)) {
    
    assert_that(ncol(y) == 1)
    y <- y[[names(y)]]
  }
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
  
  xgb <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = cv$best_iteration,
    verbose = FALSE, ...
  )
  attr(xgb, "binary") <- binary
  
  xgb
}

pred_xgb_surv <- function(xgb, df_test, intervention = NULL, X = "X") {
  
  if (!is.null(intervention)) {
    
    df_test[[X]] <- intervention
  }
  
  predict(xgb, as.matrix(df_test))
}

cross_fit_surv <- function(data, X, Z, W, time_var, event_var, time_interest, ...) {
  
  if (length(Z) == 0 & length(W) == 0) {
    
    px <- rep(mean(data[[X]]), nrow(data))
    return(list(px_z = list(1-px, px), px_zw = list(1-px, px)))
  }
  
  # determine if the setting is competing risks
  is_cr <- if (!all(data[[event_var]] %in% c(0, 1))) TRUE else FALSE
  nlvls <- max(data[[event_var]])
  
  # split into K folds
  n <- nrow(data)
  K <- 5
  folds <- sample(x = rep(1:K, each = ceiling(n / K)))[seq_len(n)]
  
  # create a list of folds to be passed to
  fld_lst <- lapply(
    seq_len(K), function(k) {
      shf <- shift(seq_len(K), n = k-1, type = "cyclic")
      tst <- folds == shf[1]
      val <- folds == shf[2]
      dev <- folds %in% shf[seq.int(3, K)]
      list(dev = dev, val = val, tst = tst)
    }
  )
  
  # outcome S/CIF via your rfs learner
  xzw_rhs <- paste(c(X, Z, W), collapse = "+")
  s_xzw_mod <- chf_rfs_cf(
    data = data[, c(X, Z, W, time_var, event_var), with = FALSE],
    X = X, time_var = time_var, event_var = event_var,
    rhs = xzw_rhs, time_interest = time_interest,
    balance_groups = FALSE, split_forest = FALSE, folds = fld_lst, ...
  )
  
  tgrid <- s_xzw_mod[[1]]$time_interest
  
  s_xz_mod <- chf_rfs_cf(
    data = data[, c(X, Z, time_var, event_var), with = FALSE],
    X = X, time_var = time_var, event_var = event_var,
    rhs = paste(c(X, Z), collapse = "+"), time_interest = tgrid,
    balance_groups = FALSE, split_forest = FALSE, folds = fld_lst, ...
  )
  
  # censoring G(t): flip event to "censoring event" (1 if censored)
  dG <- copy(data)
  dG[[event_var]] <- as.integer(dG[[event_var]] == 0L)
  g_xzw_mod <- chf_rfs_cf(
    data = dG[, c(X, Z, W, time_var, event_var), with = FALSE],
    X = X, time_var = time_var, event_var = event_var,
    rhs = xzw_rhs, time_interest = tgrid,
    balance_groups = FALSE, split_forest = FALSE, folds = fld_lst, ...
  )
  
  # time-resolved outcomes are needed
  if (!is_cr) {
    
    tres <- replicate(length(tgrid), rep(NA, n), simplify = FALSE)
    # risk indicator adjusted in 1(M > t) / G(t | x, z, w)
    ri_adj <- tres
  } else {
    
    tres <- replicate(length(tgrid), 
                      replicate(nlvls, rep(NA, n), simplify = FALSE), 
                      simplify = FALSE)
    # risk indicator adjusted in 1(M <= t, \delta = k) / G(t | x, z, w)
    ri_adj <- tres
  }
  
  y_xzw <- y_xz <- list(tres, tres)
  ey_nest <- list(list(tres, tres), list(tres, tres))
  
  # P(x | ...) elements to be filled
  px_z <- px_zw <- list(rep(NA, n), rep(NA, n))
  
  # x data
  x <- data[[X]]
  
  # if CR, compute the G(M | X, Z, W)
  if (is_cr) {
    
    gm_xzw <- rep(NA, n)
    G_m <- function(m, G, ggrid) {
      # m: length-n vector
      # G: n x K matrix, columns correspond to ggrid (length K)
      # ggrid: length-K increasing, typically starts at 0
      
      stopifnot(is.numeric(m), is.numeric(ggrid), is.matrix(G))
      stopifnot(length(ggrid) == ncol(G), length(m) == nrow(G))
      
      # clamp to grid range
      m0 <- pmax(ggrid[1], pmin(m, ggrid[length(ggrid)]))
      
      # find interval index i s.t. ggrid[i] <= m0 <= ggrid[i+1]
      i <- .bincode(m0, ggrid)              # returns in 1..length(ggrid)
      i <- pmin(i, length(ggrid) - 1L)      # cap at K-1 for i+1 indexing
      
      t0 <- ggrid[i]
      t1 <- ggrid[i + 1L]
      
      # row-wise gather of G at i and i+1
      g0 <- G[cbind(seq_along(m0), i)]
      g1 <- G[cbind(seq_along(m0), i + 1L)]
      
      # linear interpolation weight
      w <- (m0 - t0) / (t1 - t0)
      w[!is.finite(w)] <- 0                 # handles t1==t0 (shouldn't happen)
      
      (1 - w) * g0 + w * g1
    }
    
    for (i in seq_len(K)) {
      tst <- fld_lst[[i]][["tst"]]
      gm_xzw[tst] <- G_m(data[[time_var]][tst], g_xzw_mod[[i]][["srv_tst"]], tgrid)
    }
  }
  
  #' * replacing nuissance functions with GT? * 
  true_nuiss <- FALSE
  
  # cross-fit
  for (i in seq_len(K)) {
    
    # split into dev, val, tst
    tst <- fld_lst[[i]][["tst"]]
    dev <- fld_lst[[i]][["dev"]]
    val <- fld_lst[[i]][["val"]]
    
    for (t in seq_along(tgrid)) {
      
      if (!is_cr) {
        
        ri_adj[[t]][tst] <- 
          (data[[time_var]][tst] > tgrid[t]) / g_xzw_mod[[i]][["srv_tst"]][, t]
        
        for (xy in c(0, 1)) {
          
          y_xzw[[xy + 1]][[t]][tst] <- 
            s_xzw_mod[[i]][[paste0("srvx", xy, "_tst")]][, t]
          y_xz[[xy + 1]][[t]][tst] <- 
            s_xz_mod[[i]][[paste0("srvx", xy, "_tst")]][, t]
        }
      } else {
        
        for (j in seq_len(nlvls)) {
          
          ri_adj[[t]][[j]][tst] <- 
            (data[[time_var]][tst] <= tgrid[t] & data[[event_var]][tst] == j) / 
            gm_xzw[tst]
          
          for (xy in c(0, 1)) {
            
            y_xzw[[xy + 1]][[t]][[j]][tst] <- 
              s_xzw_mod[[i]][[paste0("cifx", xy, "_tst")]][, t, j]
            y_xz[[xy + 1]][[t]][[j]][tst] <- 
              s_xz_mod[[i]][[paste0("cifx", xy, "_tst")]][, t, j]
          }
        }
      }
    }
    
    # develop models on dev
    if (length(Z) > 0) {

      mod_x_z <- cv_xgb_surv(data[dev, Z, with=F], data[dev, X, with=F], ...)
    }
    
    if (length(W) > 0) {
      
      mod_x_zw <- cv_xgb_surv(data[dev, c(Z, W), with=F], data[dev, X, with=F], ...)
    } else {
      
      # inherit from Z if W empty
      mod_x_zw <- mod_x_z
    }
    
    # get the val set predictions (needed for nested means)
    px_zw_val <- pred_xgb_surv(mod_x_zw, data[val, c(Z, W), with=F])
    px_zw_val <- list(1 - px_zw_val, px_zw_val)
    
    if (length(Z) > 0) {
      
      px_z_val <- pred_xgb_surv(mod_x_z, data[val, Z, with=F])
      px_z_val <- list(1 - px_z_val, px_z_val)
    }
    
    # get the test set values
    px_zw_tst <- pred_xgb_surv(mod_x_zw, data[tst, c(Z, W), with=F])
    px_zw[[1 + 0]][tst] <- 1 - px_zw_tst
    px_zw[[1 + 1]][tst] <- px_zw_tst
    
    if (length(Z) > 0) {
      
      px_z_tst <- pred_xgb_surv(mod_x_z, data[tst, Z, with=F])
      px_z[[1 + 0]][tst] <- 1 - px_z_tst
      px_z[[1 + 1]][tst] <- px_z_tst
    } else {
      
      px_z[[1 + 0]][tst] <- 1 - mean(x)
      px_z[[1 + 1]][tst] <- mean(x)
    }
    
    if (true_nuiss) next
    for (t in seq_along(tgrid)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      
      if (!is_cr) {
        
        if (xw == xy) {
          
          ey_nest[[xw+1]][[xy+1]][[t]][tst] <- y_xz[[xy+1]][[t]][tst]
        } else {
          
          # re-fitting needed
          y_tilde <- s_xz_mod[[i]][[paste0("srvx", xy, "_val")]][, t]
          mod_nested <- cv_xgb_surv(data[val, c(X, Z), with=F], y_tilde, ...)
          ey_nest[[xw+1]][[xy+1]][[t]][tst] <- 
            pred_xgb_surv(mod_nested, data[tst, c(X, Z), with=F],
                          intervention = xw, X = X)
        }
      } else {
        
        for (j in seq_len(nlvls)) {
          
          if (xw == xy) {
            
            ey_nest[[xw+1]][[xy+1]][[t]][[j]][tst] <- y_xz[[xy+1]][[t]][[j]][tst]
          } else {
            
            # re-fitting needed
            y_tilde <- s_xz_mod[[i]][[paste0("cifx", xy, "_val")]][, t, j]
            mod_nested <- cv_xgb_surv(data[val, c(X, Z), with=F], y_tilde, ...)
            ey_nest[[xw+1]][[xy+1]][[t]][[j]][tst] <- 
              pred_xgb_surv(mod_nested, data[tst, c(X, Z), with=F],
                            intervention = xw, X = X)
          }
        }
      }
    }
  }
  
  
  if (true_nuiss) {
    
    x1_ind <- x == 1
    ey_nest_lst <- list(list(list(), list()), list(list(), list()))

    # get ground truth for S, G
    s_true <- S_T_potential_curves_from_gen(g, tgrid)
    g_true <- ground_truth_G(g, tgrid)
    
    ey_nest_lst[[0 + 1]][[0 + 1]] <- E_S_x0_given_x0_Z(g, tgrid)
    ey_nest_lst[[0 + 1]][[1 + 1]] <- E_S_x1_given_x0_Z(g, tgrid)
    ey_nest_lst[[1 + 1]][[0 + 1]] <- E_S_x0_given_x1_Z(g, tgrid)
    ey_nest_lst[[1 + 1]][[1 + 1]] <- E_S_x1_given_x1_Z(g, tgrid)
    
    for (t in seq_along(tgrid)) {
      
      ri_adj[[t]] <- (data[[time_var]] > tgrid[t]) / g_true[, t]
      
      for (xy in c(0, 1)) {
        
        y_xzw[[xy + 1]][[t]][x1_ind] <- 
          s_true[[paste0("S_x", xy, "_wx1")]][x1_ind, t]
        y_xzw[[xy + 1]][[t]][!x1_ind] <- 
          s_true[[paste0("S_x", xy, "_wx0")]][!x1_ind, t]
        
        for (xw in c(0, 1)) {
          
          ey_nest[[xw+1]][[xy+1]][[t]] <- ey_nest_lst[[xw+1]][[xy+1]][, t]
        }
      }
    }
  }
  
  list(
    ria = ri_adj,
    y_xzw = y_xzw,
    y_xz = y_xz,
    px_z = px_z,
    px_zw = px_zw,
    ey_nest = ey_nest,
    tgrid = tgrid,
    is_cr = is_cr,
    nlvls = nlvls,
    tres = tres
  )
}

pso_diff_surv <- function(cfit, data, X, Z, W, time_var, event_var, ...) {
  
  n <- nrow(data)
  
  # un-nest the cfit object
  ria <- cfit$ria
  y_xzw <- cfit$y_xzw
  y_xz <- cfit$y_xz
  px_z <- cfit$px_z
  px_zw <- cfit$px_zw
  ey_nest <- cfit$ey_nest
  tgrid <- cfit$tgrid
  is_cr <- cfit$is_cr
  nlvls <- cfit$nlvls
  
  # get data and psos
  x <- data[[X]]
  pso <- list(list(list(list(), list()), list(list(), list())),
              list(list(list(), list()), list(list(), list())))
  
  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
    
    pso[[xz+1]][[xw+1]][[xy+1]] <- cfit$tres
  }
  
  for (t in seq_along(tgrid)) {
    
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      
      # if (xy == 1 & xw == 0 & t == 20) browser()
      
      if (!is_cr) {
        
        pso[[xz+1]][[xw+1]][[xy+1]][[t]] <-
          
          # Term T1
          (x == xy) * (ria[[t]] - y_xzw[[xy+1]][[t]]) * 
          px_zw[[xw+1]] / px_zw[[xy+1]] * px_z[[xz+1]] / px_z[[xw+1]] * 
          1 / mean(x == xz) +
          
          # Term T2
          (x == xw) / mean(x == xz) * px_z[[xz+1]] / px_z[[xw+1]] *
          (y_xzw[[xy+1]][[t]] - ey_nest[[xw+1]][[xy+1]][[t]]) +
          
          # Term T3
          (x == xz) / mean(x == xz) * ey_nest[[xw+1]][[xy+1]][[t]]
      } else {
        
        for (j in seq_len(nlvls)) {
          
          pso[[xz+1]][[xw+1]][[xy+1]][[t]][[j]] <-
            
            # Term T1
            (x == xy) * (ria[[t]][[j]] - y_xzw[[xy+1]][[t]][[j]]) * 
            px_zw[[xw+1]] / px_zw[[xy+1]] * px_z[[xz+1]] / px_z[[xw+1]] * 
            1 / mean(x == xz) +
            
            # Term T2
            (x == xw) / mean(x == xz) * px_z[[xz+1]] / px_z[[xw+1]] *
            (y_xzw[[xy+1]][[t]][[j]] - ey_nest[[xw+1]][[xy+1]][[t]][[j]]) +
            
            # Term T3
            (x == xz) / mean(x == xz) * ey_nest[[xw+1]][[xy+1]][[t]][[j]]
        }
      }
    }
  }
  
  pso
}

measure_spec <- function() {
  
  list(
    tv = list(
      sgn = c(1, -1),
      spc = list(c(1, 1, 1), c(0, 0, 0)),
      nm = "tv"
    ),
    ctfde = list(
      sgn = c(1, -1),
      spc = list(c(0, 0, 1), c(0, 0, 0)),
      nm = "ctfde"
    ),
    ctfie = list(
      sgn = c(1, -1),
      spc = list(c(0, 0, 1), c(0, 1, 1)),
      nm = "ctfie"
    ),
    ctfse = list(
      sgn = c(1, -1),
      spc = list(c(0, 1, 1), c(1, 1, 1)),
      nm = "ctfse"
    )
  )
}

one_step_debias_surv <- function(data, X, Z, W, time_var, event_var, 
                                 time_interest = NULL, eps_trim = 0, 
                                 copula = NULL, tau_grid = 0, ...) {
  
  cfit <- cross_fit_surv(data, X, Z, W, time_var, event_var, time_interest, ...)
  pso <- pso_diff_surv(cfit, data, X, Z, W, time_var, event_var, ...)
  
  # get extreme propensity weights
  extrm_pxz <- (cfit$px_z[[1]] < eps_trim) | (1 - cfit$px_z[[1]] < eps_trim)
  extrm_pxzw <-  (cfit$px_zw[[1]] < eps_trim) | (1 - cfit$px_zw[[1]] < eps_trim)
  extrm_idx <- extrm_pxz | extrm_pxzw
  
  # report if large number of propensity weights below specified threshold
  if (mean(extrm_idx) > 0.02) {
    message(round(100 * mean(extrm_idx), 2),
            "% of extreme P(x | z) or P(x | z, w) probabilities at threshold",
            " = ", eps_trim, ".\n",
            "Reported results are for the overlap population. ",
            "Consider investigating overlap issues.")
  }
  
  # trim population to extreme weights
  if (!cfit$is_cr) {
    
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1))
      for (t in seq_along(cfit$tgrid))
        pso[[xz+1]][[xw+1]][[xy+1]][[t]][extrm_idx] <- NA
  } else {
    
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) 
      for (t in seq_along(cfit$tgrid))
        for (j in seq_len(cfit$nlvls))
          pso[[xz+1]][[xw+1]][[xy+1]][[t]][[j]][extrm_idx] <- NA
  }
  
  # get specification of measures to be reported
  eff <- measure_spec()
  
  is_sens <- if (!is.null(copula)) TRUE else FALSE
  if (is_sens) {
    
    trim_cif <- function(cif) {
      
      cif <- pmin(pmax(cif, 0), 1)
      cummax(cif)
    }

    res <- c()
    elm <- list(
      list(mu = rep(NA, length(cfit$tgrid)), sd = rep(NA, length(cfit$tgrid))),
      list(mu = rep(NA, length(cfit$tgrid)), sd = rep(NA, length(cfit$tgrid)))
    )
    cif <- list(list(list(elm, elm), list(elm, elm)), 
                list(list(elm, elm), list(elm, elm))) 
    
    elm2 <- replicate(length(tau_grid), list(mean = rep(NA, length(cfit$tgrid)),
                              lwr = rep(NA, length(cfit$tgrid)),
                              upr = rep(NA, length(cfit$tgrid))), 
                      simplify = FALSE)
    shat <- list(list(list(elm2, elm2), list(elm2, elm2)), 
                 list(list(elm2, elm2), list(elm2, elm2))) 
    
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1))
      for (t in seq_along(cfit$tgrid))
        for (j in seq_len(cfit$nlvls)) {
          
          pseudo_out <- pso[[xz+1]][[xw+1]][[xy+1]][[t]][[j]]
          psi_osd <- mean(pseudo_out, na.rm = TRUE)
          dev <- sqrt(var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))
          
          
          cif[[xz+1]][[xw+1]][[xy+1]][[j]][["mu"]][t] <- psi_osd 
          cif[[xz+1]][[xw+1]][[xy+1]][[j]][["sd"]][t] <- dev
        }
    
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      
      cif_a <- cif[[xz+1]][[xw+1]][[xy+1]][[1]][["mu"]]
      sd_a <- cif[[xz+1]][[xw+1]][[xy+1]][[1]][["sd"]]
      
      cif_b <- cif[[xz+1]][[xw+1]][[xy+1]][[2]][["mu"]]
      sd_b <- cif[[xz+1]][[xw+1]][[xy+1]][[2]][["sd"]]
      
      for (tau_id in seq_along(tau_grid)) {
        
        shat_mean <- cif_copula_single(cif_a, cif_b, copula, tau_grid[tau_id])
        shat[[xz+1]][[xw+1]][[xy+1]][[tau_id]][["mean"]] <- shat_mean
        
        shat_upr <- cif_copula_single(trim_cif(cif_a - 1.96 * sd_a), 
                                      trim_cif(cif_b + 1.96 * sd_b), 
                                      copula, tau_grid[tau_id])
        shat[[xz+1]][[xw+1]][[xy+1]][[tau_id]][["upr"]] <- shat_upr
        
        shat_lwr <- cif_copula_single(trim_cif(cif_a + 1.96 * sd_a), 
                                      trim_cif(cif_b - 1.96 * sd_b), 
                                      copula, tau_grid[tau_id])
        if (tau_grid[tau_id] == 0.5) browser()
        shat[[xz+1]][[xw+1]][[xy+1]][[tau_id]][["lwr"]] <- shat_lwr
      }
    }
    
    for (i in seq_along(eff)) {
      
      xz1 <- eff[[i]]$spc[[1]][1]
      xw1 <- eff[[i]]$spc[[1]][2]
      xy1 <- eff[[i]]$spc[[1]][3]
      
      xz2 <- eff[[i]]$spc[[2]][1]
      xw2 <- eff[[i]]$spc[[2]][2]
      xy2 <- eff[[i]]$spc[[2]][3]
      
      # get the point estimates for tau
      for (tau_id in seq_along(tau_grid)) {
        
        eff_mean <- shat[[xz1+1]][[xw1+1]][[xy1+1]][[tau_id]][["mean"]] -
                    shat[[xz2+1]][[xw2+1]][[xy2+1]][[tau_id]][["mean"]]
        
        eff_lwr <- shat[[xz1+1]][[xw1+1]][[xy1+1]][[tau_id]][["lwr"]] -
                   shat[[xz2+1]][[xw2+1]][[xy2+1]][[tau_id]][["upr"]]
        
        eff_upr <- shat[[xz1+1]][[xw1+1]][[xy1+1]][[tau_id]][["upr"]] -
                   shat[[xz2+1]][[xw2+1]][[xy2+1]][[tau_id]][["lwr"]]
        
        res <- rbind(
          res,
          data.frame(effect = eff[[i]]$nm, value = eff_mean,
                     lwr = eff_lwr, upr = eff_upr, tau = tau_grid[tau_id],
                     time_interest = cfit$tgrid)
        )
      }
    }
    
    return(res)
  }
  
  res <- c()
  for (i in seq_along(eff)) {
    
    for (t in seq_along(cfit$tgrid)) {
      
      if (!cfit$is_cr) {
        
        pseudo_out <- 0
        for (s in seq_along(eff[[i]]$sgn)) {
          
          xz <- eff[[i]]$spc[[s]][1]
          xw <- eff[[i]]$spc[[s]][2]
          xy <- eff[[i]]$spc[[s]][3]
          pseudo_out <- pseudo_out + eff[[i]]$sgn[s] * pso[[xz+1]][[xw+1]][[xy+1]][[t]]
        }
        psi_osd <- mean(pseudo_out, na.rm = TRUE)
        dev <- sqrt(var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))
        
        res <- rbind(
          res,
          data.frame(effect = eff[[i]]$nm, value = psi_osd, sd = dev, 
                     time_interest = cfit$tgrid[t])
        ) 
      } else {
        
        for (j in seq_len(cfit$nlvls)) {
          
          pseudo_out <- 0
          for (s in seq_along(eff[[i]]$sgn)) {
            
            xz <- eff[[i]]$spc[[s]][1]
            xw <- eff[[i]]$spc[[s]][2]
            xy <- eff[[i]]$spc[[s]][3]
            pseudo_out <- pseudo_out + 
              eff[[i]]$sgn[s] * pso[[xz+1]][[xw+1]][[xy+1]][[t]][[j]]
          }
          psi_osd <- mean(pseudo_out, na.rm = TRUE)
          dev <- sqrt(var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))
          
          res <- rbind(
            res,
            data.frame(effect = eff[[i]]$nm, value = psi_osd, sd = dev, 
                       time_interest = cfit$tgrid[t], event = j)
          )  
        }
      }
    }
  }
  
  # pw <- list(px_z = cfit$px_z, px_zw = cfit$px_zw)
  # attr(res, "pw") <- pw
  as.data.table(res)
}
