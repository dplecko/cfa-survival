
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
    nrounds = cv$early_stop$best_iteration,
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

cross_fit_surv <- function(data, X, Z, W, time_var, event_var, time_interest, 
                           martingale_debias, modify_S, modify_G, dgm, ...) {
  
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
  if (is.null(modify_S)) {
    s_xzw_mod <- chf_rfs_cf(
      data = data[, c(X, Z, W, time_var, event_var), with = FALSE],
      X = X, time_var = time_var, event_var = event_var,
      rhs = xzw_rhs, time_interest = time_interest,
      balance_groups = FALSE, split_forest = FALSE, folds = fld_lst, ...
    )
    tgrid <- s_xzw_mod[[1]]$time_interest
  } else {
    
    assert_that(!is.null(time_interest), 
                msg = "Oracle runs must have `time_interest` specified.")
    tgrid <- time_interest
    s_xzw_mod <- vector("list", K)
  }
  
  # censoring G(t): flip event to "censoring event" (1 if censored)
  dG <- copy(data)
  dG[[event_var]] <- as.integer(dG[[event_var]] == 0L)
  
  if (is.null(modify_G)) {
    
    g_xzw_mod <- chf_rfs_cf(
      data = dG[, c(X, Z, W, time_var, event_var), with = FALSE],
      X = X, time_var = time_var, event_var = event_var,
      rhs = xzw_rhs, time_interest = tgrid,
      balance_groups = FALSE, split_forest = FALSE, folds = fld_lst, ...
    )
  } else {
    
    g_xzw_mod <- vector("list", K)
  }
  
  # ---- nuisance modification: marginal KM corruption or oracle truth ----
  # dgm is only ever available for synthetic experiments (known DGP); on
  # real data (no dgm) neither modify_S/modify_G nor the A0/A1/A2 oracle
  # decomposition are meaningful, so this whole block is skipped there.
  ora_px_z <- ora_px_zw <- ora_nu <- NULL
  if ((!is.null(modify_S) || !is.null(modify_G)) && is.null(dgm)) {
    stop("modify_S/modify_G require `dgm` (the known data-generating ",
         "mechanism) to compute oracle/corrupt nuisances; got dgm = NULL.")
  }
  if (!is.null(dgm)) {

    # first computing the corruptions / ground truth
    if (!is.null(modify_S)) assert_that(modify_S %in% c("corrupt", "oracle"))
    if (!is.null(modify_G)) assert_that(modify_G %in% c("corrupt", "oracle"))

    km_eval <- function(times, events, tgrid) {
      sf <- survival::survfit(survival::Surv(times, events) ~ 1)
      fn <- stepfun(sf$time, c(1, sf$surv), right = FALSE)
      pmax(fn(tgrid), 1e-6)
    }

    # oracle S(t | X,Z,W) / CIF_k(t | X,Z,W), and their counterfactuals
    # under X:=x (unit's own factual Z,W throughout, only X intervened)
    if (is_cr) {
      
      x1_ind <- data[[X]] == 1
      cifx0 <- array(0, dim = c(n, length(tgrid), nlvls))
      cifx1 <- array(0, dim = c(n, length(tgrid), nlvls))
      
      for (j in seq_len(nlvls)) {
        c_true <- CIF_conditional_exact(dgm, tgrid, event = j)
        cifx0[, , j] <- c_true$cifx0
        cifx1[, , j] <- c_true$cifx1
      }
      
      cif <- cifx0
      cif[x1_ind, , ] <- cifx1[x1_ind, , ]
      
      srvx0 <- 1 - apply(cifx0, c(1, 2), sum)
      srvx1 <- 1 - apply(cifx1, c(1, 2), sum)
      srv <- srvx0
      srv[x1_ind, ] <- srvx1[x1_ind, ]
      
      # closed-form oracle nuisances for the A0/A1/A2 decomposition
      ora_px_z  <- oracle_px_z(dgm)
      ora_px_zw <- oracle_px_zw(dgm, ora_px_z)
      ora_nu    <- oracle_nuisances_cr(dgm, tgrid)$ey_nest

    } else {
      
      x1_ind <- data[[X]] == 1
      s_true <- S_T_potential_curves_from_gen(dgm, tgrid)
      
      srvx0 <- s_true[["S_x0_wx0"]]
      srvx0[x1_ind, ] <- s_true[["S_x0_wx1"]][x1_ind, ]
      
      srvx1 <- s_true[["S_x1_wx0"]]
      srvx1[x1_ind, ] <- s_true[["S_x1_wx1"]][x1_ind, ]
      
      srv <- srvx0
      srv[x1_ind, ] <- srvx1[x1_ind, ]
    }
    
    # oracle G(t | X,Z,W)
    g_true <- ground_truth_G(dgm, tgrid)
    h_true <- -log(pmax(g_true, 1e-6))
    
    for (i in seq_len(K)) {
      
      dev <- fld_lst[[i]]$dev
      val <- fld_lst[[i]]$val
      tst <- fld_lst[[i]]$tst
      
      n_tst <- sum(tst)
      n_val <- sum(val)
      
      if (identical(modify_S, "corrupt")) {
        
        s_vec <- km_eval(data[[time_var]][dev], data[[event_var]][dev], tgrid)
        Mt <- matrix(s_vec, n_tst, length(tgrid), byrow = TRUE)
        Mv <- matrix(s_vec, n_val, length(tgrid), byrow = TRUE)
        
        s_xzw_mod[[i]]$srv_tst <- Mt
        s_xzw_mod[[i]]$srvx0_tst <- Mt
        s_xzw_mod[[i]]$srvx1_tst <- Mt
        s_xzw_mod[[i]]$srvx0_val <- Mv
        s_xzw_mod[[i]]$srvx1_val <- Mv
        
      } else if (identical(modify_S, "oracle") && is_cr) {
        
        s_xzw_mod[[i]]$srv_tst <- srv[tst, ]
        s_xzw_mod[[i]]$cif_tst <- cif[tst, , ]
        s_xzw_mod[[i]]$cifx0_tst <- cifx0[tst, , ]
        s_xzw_mod[[i]]$cifx1_tst <- cifx1[tst, , ]
        s_xzw_mod[[i]]$cifx0_val <- cifx0[val, , ]
        s_xzw_mod[[i]]$cifx1_val <- cifx1[val, , ]
        
      } else if (identical(modify_S, "oracle")) {
        
        s_xzw_mod[[i]]$srv_tst <- srv[tst, ]
        s_xzw_mod[[i]]$srvx0_tst <- srvx0[tst, ]
        s_xzw_mod[[i]]$srvx1_tst <- srvx1[tst, ]
        s_xzw_mod[[i]]$srvx0_val <- srvx0[val, ]
        s_xzw_mod[[i]]$srvx1_val <- srvx1[val, ]
      }
      
      if (identical(modify_G, "corrupt")) {
        
        g_vec <- km_eval(data[[time_var]][dev], 1L - data[[event_var]][dev], tgrid)
        Mt <- matrix(g_vec, n_tst, length(tgrid), byrow = TRUE)
        Ht <- matrix(-log(g_vec), n_tst, length(tgrid), byrow = TRUE)
        
        g_xzw_mod[[i]]$srv_tst <- Mt
        g_xzw_mod[[i]]$chf_tst <- Ht
        
      } else if (identical(modify_G, "oracle")) {
        
        g_xzw_mod[[i]]$srv_tst <- g_true[tst, ]
        g_xzw_mod[[i]]$chf_tst <- h_true[tst, ]
      }
    }
  }
  
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
  
  y_xzw <- y_xzw_ora <- list(tres, tres)
  ey_nest <- ey_nest_ora <- list(list(tres, tres), list(tres, tres))
  
  # P(x | ...) elements to be filled
  px_z <- px_zw <- list(rep(NA, n), rep(NA, n))
  
  # x data
  x <- data[[X]]
  
  # if CR, compute the G(M | X, Z, W)
  if (is_cr) {
    
    gm_xzw <- rep(NA, n)
    G_m <- function(m, G, ggrid) {
      
      stopifnot(is.numeric(m), is.numeric(ggrid), is.matrix(G))
      stopifnot(length(ggrid) == ncol(G), length(m) == nrow(G))
      
      # clamp to grid range
      m0 <- pmax(ggrid[1], pmin(m, ggrid[length(ggrid)]))
      
      # find interval index i s.t. ggrid[i] <= m0 <= ggrid[i+1]
      i <- .bincode(m0, ggrid)              
      i <- pmin(i, length(ggrid) - 1L)
      
      t0 <- ggrid[i]
      t1 <- ggrid[i + 1L]
      
      # row-wise gather of G at i and i+1
      g0 <- G[cbind(seq_along(m0), i)]
      g1 <- G[cbind(seq_along(m0), i + 1L)]
      
      # linear interpolation weight
      w <- (m0 - t0) / (t1 - t0)
      w[!is.finite(w)] <- 0
      
      (1 - w) * g0 + w * g1
    }
    
    for (i in seq_len(K)) {
      tst <- fld_lst[[i]][["tst"]]
      gm_xzw[tst] <- G_m(data[[time_var]][tst], g_xzw_mod[[i]][["srv_tst"]], tgrid)
    }
  }
  
  # cross-fit
  for (i in seq_len(K)) {
    
    # split into dev, val, tst
    tst <- fld_lst[[i]][["tst"]]
    dev <- fld_lst[[i]][["dev"]]
    val <- fld_lst[[i]][["val"]]
    
    M_tst <- data[[time_var]][tst]
    Mt <- pmax(1L, findInterval(M_tst, tgrid))
    S_M <- s_xzw_mod[[i]][["srv_tst"]][cbind(seq_along(Mt), Mt)]
    G_M <- g_xzw_mod[[i]][["srv_tst"]][cbind(seq_along(Mt), Mt)]
    xi2_int <- 0
    
    zeta2_A <- rep(0, sum(tst))
    zeta2_B <- matrix(0, nrow = sum(tst), ncol = nlvls)
    
    if (is_cr) {
      
      cif_tst <- s_xzw_mod[[i]][["cif_tst"]]
      s_all_tst <- s_xzw_mod[[i]][["srv_tst"]]
      g_tst <- g_xzw_mod[[i]][["srv_tst"]]
      h_c_tst <- g_xzw_mod[[i]][["chf_tst"]]
      
      # Grid approximation to F_j(M-) and S_all(M-)
      m_idx <- findInterval(M_tst, tgrid)
      has_m_idx <- m_idx > 0
      
      s_m_minus <- rep(1, length(M_tst))
      s_m_minus[has_m_idx] <-
        s_all_tst[cbind(which(has_m_idx), m_idx[has_m_idx])]
      
      cif_m_minus <- matrix(
        0,
        nrow = length(M_tst),
        ncol = nlvls
      )
      
      for (j in seq_len(nlvls)) {
        cif_m_minus[has_m_idx, j] <-
          cif_tst[cbind(which(has_m_idx), m_idx[has_m_idx], j)]
      }
    }
    
    for (t in seq_along(tgrid)) {
      
      if (!is_cr) {
        
        ri_adj[[t]][tst] <- 
          (data[[time_var]][tst] > tgrid[t]) / g_xzw_mod[[i]][["srv_tst"]][, t]
        
        if (martingale_debias) {
          
          S_t <- s_xzw_mod[[i]][["srv_tst"]][, t]
          G_t <- g_xzw_mod[[i]][["srv_tst"]][, t]
          H_tst <- g_xzw_mod[[i]][["chf_tst"]]
          dH <- H_tst[, t] - if (t == 1) 0 else H_tst[, t - 1]
          
          # Accumulate Term \xi_2 (Continuous penalty integral up to t)
          xi2_int <- xi2_int + (M_tst >= tgrid[t]) * dH / (S_t * G_t)
          
          # Compute Term \xi_1 (Censoring event jump)
          xi1 <- S_t * (M_tst <= tgrid[t] & data[[event_var]][tst] == 0) / (S_M * G_M)
          
          # Combine: IPCW + \xi_1(t) - \xi_2(t)
          ri_adj[[t]][tst] <- ri_adj[[t]][tst] + xi1 - (S_t * xi2_int)
        }
        
        for (xy in c(0, 1)) {
          
          y_xzw[[xy + 1]][[t]][tst] <-
            s_xzw_mod[[i]][[paste0("srvx", xy, "_tst")]][, t]
        }
      } else {
        
        if (martingale_debias) {
          
          G_t <- pmax(g_tst[, t], 1e-6)
          dH <- h_c_tst[, t] -
            if (t == 1) 0 else h_c_tst[, t - 1]
          
          S_left <- if (t == 1) {
            rep(1, sum(tst))
          } else {
            s_all_tst[, t - 1]
          }
          S_left <- pmax(S_left, 1e-6)
          
          zeta_weight <-
            (M_tst >= tgrid[t]) * dH / (G_t * S_left)
          
          zeta2_A <- zeta2_A + zeta_weight
          
          for (j in seq_len(nlvls)) {
            
            F_left <- if (t == 1) {
              rep(0, sum(tst))
            } else {
              cif_tst[, t - 1, j]
            }
            
            zeta2_B[, j] <-
              zeta2_B[, j] + zeta_weight * F_left
          }
        }
        
        for (j in seq_len(nlvls)) {
          
          F_t <- cif_tst[, t, j]
          
          # IPCW event term
          ri_adj[[t]][[j]][tst] <-
            (M_tst <= tgrid[t] &
               data[[event_var]][tst] == j) /
            pmax(gm_xzw[tst], 1e-6)
          
          if (martingale_debias) {
            
            # Censoring-event jump term
            zeta1 <-
              (M_tst <= tgrid[t] &
                 data[[event_var]][tst] == 0) /
              pmax(gm_xzw[tst], 1e-6) *
              (F_t - cif_m_minus[, j]) /
              pmax(s_m_minus, 1e-6)
            
            # Continuous censoring augmentation
            zeta2 <- F_t * zeta2_A - zeta2_B[, j]
            
            ri_adj[[t]][[j]][tst] <-
              ri_adj[[t]][[j]][tst] + zeta1 - zeta2
          }
          
          for (xy in c(0, 1)) {
            
            y_xzw[[xy + 1]][[t]][[j]][tst] <-
              s_xzw_mod[[i]][[paste0("cifx", xy, "_tst")]][, t, j]

            if (!is.null(dgm)) {
              cif_ora <- if (xy == 1) cifx1[tst, , ] else cifx0[tst, , ]
              y_xzw_ora[[xy + 1]][[t]][[j]][tst] <- cif_ora[, t, j]
            }
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
    
    # if (true_nuiss) next
    for (t in seq_along(tgrid)) for (xy in c(0, 1)) {
      
      if (!is_cr) {
        
        y_tilde <- s_xzw_mod[[i]][[paste0("srvx", xy, "_val")]][, t]
        mod_nested <- cv_xgb_surv(data[val, c(X, Z), with=F], y_tilde, ...)
        for (xw in c(0, 1))
          ey_nest[[xw+1]][[xy+1]][[t]][tst] <-
          pred_xgb_surv(mod_nested, data[tst, c(X, Z), with=F], intervention=xw, X=X)
      } else {
        
        for (j in seq_len(nlvls)) {
          
          
          # oracle ey_nest (A0/A1/A2 decomposition): skipped when modify_S is
          # already "oracle", since ey_nest itself is then already a fit on
          # the oracle target and a separate ey_nest_ora would just be a
          # redundant second fit of the identical (data, target) pair
          if (!is.null(dgm) && !identical(modify_S, "oracle")) {

            cif_ora <- if (xy == 1) cifx1[val, ,] else cifx0[val, ,]
            y_tilde_ora <- cif_ora[, t, j]
            mod_nested_ora <- cv_xgb_surv(data[val, c(X, Z), with=F], y_tilde_ora, ...)
            for (xw in c(0, 1))
              ey_nest_ora[[xw+1]][[xy+1]][[t]][[j]][tst] <-
              pred_xgb_surv(mod_nested_ora, data[tst, c(X, Z), with=F], intervention=xw, X=X) 
          }
          
          # estimated ey_nest
          y_tilde <- s_xzw_mod[[i]][[paste0("cifx", xy, "_val")]][, t, j]
          mod_nested <- cv_xgb_surv(data[val, c(X, Z), with=F], y_tilde, ...)
          for (xw in c(0, 1))
            ey_nest[[xw+1]][[xy+1]][[t]][[j]][tst] <-
            pred_xgb_surv(mod_nested, data[tst, c(X, Z), with=F], intervention=xw, X=X)
        }
      }
    }
    
  }
  
  list(
    ria = ri_adj,
    y_xzw = y_xzw,
    y_xzw_ora = y_xzw_ora,
    px_z = px_z,
    px_zw = px_zw,
    ey_nest = ey_nest,
    ey_nest_ora = ey_nest_ora,
    ora_px_z = ora_px_z,
    ora_px_zw = ora_px_zw,
    ora_nu = ora_nu,
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
  y_xzw_ora <- cfit$y_xzw_ora
  px_z <- cfit$px_z
  px_zw <- cfit$px_zw
  ey_nest <- cfit$ey_nest
  ey_nest_ora <- cfit$ey_nest_ora
  tgrid <- cfit$tgrid
  is_cr <- cfit$is_cr
  nlvls <- cfit$nlvls
  
  # get data and psos
  x <- data[[X]]
  pso <- list(list(list(list(), list()), list(list(), list())),
              list(list(list(), list()), list(list(), list())))
  pso_ora <- pso
  
  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
    
    pso[[xz+1]][[xw+1]][[xy+1]] <- cfit$tres
    pso_ora[[xz+1]][[xw+1]][[xy+1]] <- cfit$tres
  }
  
  for (t in seq_along(tgrid)) {
    
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {

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
          
          # oracle pseudo-outcome phi(V; P) + psi * 1(x = xz) / P_hat(xz):
          # same ria (exact under modify_G = "oracle", martingale off),
          # oracle CIF, oracle propensities, closed-form oracle nu
          if (!is.null(cfit$ora_px_z)) {
            
            oz <- cfit$ora_px_z; ozw <- cfit$ora_px_zw
            onu <- cfit$ora_nu[[xw+1]][[xy+1]][[t]][[j]]
            
            pso_ora[[xz+1]][[xw+1]][[xy+1]][[t]][[j]] <-
              (x == xy) * (ria[[t]][[j]] - y_xzw_ora[[xy+1]][[t]][[j]]) *
              ozw[[xw+1]] / ozw[[xy+1]] * oz[[xz+1]] / oz[[xw+1]] *
              1 / mean(x == xz) +
              (x == xw) / mean(x == xz) * oz[[xz+1]] / oz[[xw+1]] *
              (y_xzw_ora[[xy+1]][[t]][[j]] - onu) +
              (x == xz) / mean(x == xz) * onu
          }
        }
      }
    }
  }
  
  attr(pso, "pso_ora") <- pso_ora
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
                                 copula = NULL, tau_grid = 0, 
                                 martingale_debias = TRUE, 
                                 modify_S = NULL, modify_G = NULL, 
                                 dgm = NULL, ...) {
  
  cfit <- cross_fit_surv(data, X, Z, W, time_var, event_var, time_interest, 
                         martingale_debias, modify_S, modify_G, dgm, ...)
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
  res_sens <- NULL
  if (is_sens) {
    
    res_sens <- c()
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
      
      # four corners + sampling
      for (tau_id in seq_along(tau_grid)) {
        
        sh <- cuatro_esquinas(
          c1 = cif_a, s1 = sd_a,
          c2 = cif_b, s2 = sd_b,
          copula = copula, tau = tau_grid[tau_id],
          z = 1.96, trim_cif = trim_cif,
          do_samp = TRUE
        )
        
        shat[[xz+1]][[xw+1]][[xy+1]][[tau_id]][["mean"]] <- sh[["mean"]]
        shat[[xz+1]][[xw+1]][[xy+1]][[tau_id]][["lwr"]]  <- sh[["lwr"]]
        shat[[xz+1]][[xw+1]][[xy+1]][[tau_id]][["upr"]]  <- sh[["upr"]]
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
        
        res_sens <- rbind(
          res_sens,
          data.frame(effect = eff[[i]]$nm, value = eff_mean,
                     lwr = eff_lwr, upr = eff_upr, tau = tau_grid[tau_id],
                     time_interest = cfit$tgrid)
        )
      }
    }
  }
  
  { # classical point estimates (always computed, cheap relative to the fit)
    
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
  }
  
  { # raw potential-outcome cells psi(xz,xw,xy), pre-differencing -- just a
    # summary of pso (already computed above), useful for isolating which
    # nested counterfactual a bias comes from without re-deriving effects
    res_po <- c()
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      
      for (t in seq_along(cfit$tgrid)) {
        
        if (!cfit$is_cr) {
          
          pseudo_out <- pso[[xz+1]][[xw+1]][[xy+1]][[t]]
          psi_osd <- mean(pseudo_out, na.rm = TRUE)
          dev <- sqrt(var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))
          
          res_po <- rbind(
            res_po,
            data.frame(xz = xz, xw = xw, xy = xy, value = psi_osd, sd = dev,
                       time_interest = cfit$tgrid[t])
          )
        } else {
          
          for (j in seq_len(cfit$nlvls)) {
            
            pseudo_out <- pso[[xz+1]][[xw+1]][[xy+1]][[t]][[j]]
            psi_osd <- mean(pseudo_out, na.rm = TRUE)
            dev <- sqrt(var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))
            
            # ---- A0/A1/A2 decomposition columns ----
            # A0 = value_ora - truth (merge with ground_truth_cr_cells downstream)
            # A2 ~= e2_t2 + e2_t3 (empirical E2, censoring term = 0 under oracle G)
            # A1 = (value - value_ora) - A2
            value_ora <- e2_t2 <- e2_t3 <- a1 <- NA_real_
            pso_ora <- attr(pso, "pso_ora")
            if (!is.null(cfit$ora_px_z)) {
              
              xvec <- data[[X]]
              value_ora <- mean(pso_ora[[xz+1]][[xw+1]][[xy+1]][[t]][[j]],
                                na.rm = TRUE)
              
              oz <- cfit$ora_px_z; ozw <- cfit$ora_px_zw
              err_y <- cfit$y_xzw_ora[[xy+1]][[t]][[j]] -
                cfit$y_xzw[[xy+1]][[t]][[j]]
              err_nu <- cfit$ora_nu[[xw+1]][[xy+1]][[t]][[j]] -
                cfit$ey_nest[[xw+1]][[xy+1]][[t]][[j]]
              
              # r1 = lam_hat / lam on Z-scale; r2 = lam / lam_hat on ZW-scale
              r1 <- (cfit$px_z[[xz+1]] / cfit$px_z[[xw+1]]) /
                (oz[[xz+1]] / oz[[xw+1]])
              r2 <- (ozw[[xy+1]] / ozw[[xw+1]]) /
                (cfit$px_zw[[xy+1]] / cfit$px_zw[[xw+1]])
              wzw <- (xvec == xw) * oz[[xz+1]] / oz[[xw+1]] / mean(xvec == xz)
              
              e2_t2 <- mean(wzw * r1 * (r2 - 1) * err_y, na.rm = TRUE)
              e2_t3 <- mean((xvec == xz) / mean(xvec == xz) * (r1 - 1) * err_nu,
                            na.rm = TRUE)
              a1 <- (psi_osd - value_ora) - e2_t2 - e2_t3
            }
            
            res_po <- rbind(
              res_po,
              data.frame(xz = xz, xw = xw, xy = xy, value = psi_osd, sd = dev,
                         value_ora = value_ora, a1 = a1,
                         e2_t2 = e2_t2, e2_t3 = e2_t3,
                         time_interest = cfit$tgrid[t], event = j)
            )
          }
        }
      }
    }
  }
  
  structure(
    list(
      measures = as.data.table(res),
      measures_sens = if (!is.null(res_sens)) as.data.table(res_sens) else NULL,
      measures_po = as.data.table(res_po),
      is_cr = cfit$is_cr, is_sens = is_sens,
      time_interest = cfit$tgrid, copula = copula, tau_grid = tau_grid
    ), class = "fairsurv_osd"
  )
}

autoplot.fairsurv_osd <- function(object, ...) {
  
  meas <- c("ctfde", "ctfie", "ctfse", "tv")
  
  alpha <- 0.05 # fixed for now
  width <- qnorm(1 - alpha / 2)
  
  plt_dat <- copy(if (object$is_sens) object$measures_sens else object$measures)
  plt_dat[, effect := factor(effect, levels = c("tv", "ctfde", "ctfie", "ctfse"),
                             labels = c("Total Variation", "Direct", "Indirect", "Spurious"))]
  if (object$is_sens) { #
    
    # tau as linetype (assumes up to 4 tau values)
    plt_dat[, tau_f := factor(tau, levels = sort(unique(tau)))]
    lt_vals <- c("solid", "dashed", "dotted", "dotdash")[seq_len(length(levels(plt_dat$tau_f)))]
    
    # union envelope over tau (one ribbon per effect)
    env <- plt_dat[, .(lwr = min(lwr), upr = max(upr)),
                   by = c("effect", "time_interest")]
    
    p <- ggplot(plt_dat, aes(x = time_interest, y = value, color = effect)) +
      geom_ribbon(
        data = env,
        aes(ymin = lwr, ymax = upr, fill = effect, x = time_interest),
        alpha = 0.2, linewidth = 0, inherit.aes = FALSE
      ) +
      geom_line(aes(linetype = tau_f), linewidth = 1) +
      theme_bw(base_size = 14) +
      xlab("Time") + ylab("Effect Value") +
      scale_color_discrete(
        name = "Effect", labels = c("TV", "x-DE", "x-IE", "x-SE")
      ) +
      scale_fill_discrete(
        name = "Effect", labels = c("TV", "x-DE", "x-IE", "x-SE")) +
      scale_linetype_manual(name = latex2exp::TeX("$\\tau$"), values = lt_vals) +
      facet_wrap(~effect, ncol = 4, scales = "free") +
      theme(legend.position = "bottom", axis.text = element_text(size = 9)) + 
      guides(
        linetype = guide_legend(
          title.theme = element_text(size = 20)
        )
      ) +
      scale_y_continuous(labels = scales::percent)
    
  } else if (object$is_cr) { # competing risks
    
    p <- ggplot(plt_dat, aes(x = time_interest, y = value, color = effect, 
                             fill = effect)) +
      geom_line() + theme_bw() +
      geom_ribbon(aes(ymin = value - width * sd, ymax = value + width * sd),
                  alpha = 0.4, linewidth = 0) +
      facet_wrap(~event, scales = "free")
  } else {
    
    p <- ggplot(plt_dat, aes(x = time_interest, y = value, color = effect, 
                             fill = effect)) +
      geom_line(linewidth = 1) + theme_bw(base_size = 14) +
      geom_ribbon(aes(ymin = value - width * sd, ymax = value + width * sd),
                  alpha = 0.3, linewidth = 0) +
      facet_wrap(~effect, ncol = 4, scales = "free") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_discrete(
        name = "Effect",
        labels = c("TV", "x-DE", "x-IE", "x-SE")
      ) +
      scale_color_discrete(
        name = "Effect",
        labels = c("TV", "x-DE", "x-IE", "x-SE")
      ) +
      xlab("Time (days)") + ylab("Effect Value") +
      theme(legend.position = "bottom", axis.text = element_text(size = 9))
  }
  
  p
}

cuatro_esquinas <- function(c1, s1, c2, s2, copula, tau, z = 1.96, trim_cif,
                            do_samp = FALSE, B = 200, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  c1_l <- trim_cif(c1 - z * s1); c1_u <- trim_cif(c1 + z * s1)
  c2_l <- trim_cif(c2 - z * s2); c2_u <- trim_cif(c2 + z * s2)
  
  # keep CIF_1 + CIF_2 <= 1 (cheap projection)
  proj_sum <- function(a, b) {
    b <- pmin(b, 1 - a)
    a <- pmin(a, 1 - b)
    list(a = trim_cif(a), b = trim_cif(b))
  }
  
  sh_mu <- cif_copula_single(c1, c2, copula, tau)
  
  crn <- list(
    proj_sum(c1_l, c2_l),
    proj_sum(c1_l, c2_u),
    proj_sum(c1_u, c2_l),
    proj_sum(c1_u, c2_u)
  )
  sh1 <- cif_copula_single(crn[[1]]$a, crn[[1]]$b, copula, tau)
  sh2 <- cif_copula_single(crn[[2]]$a, crn[[2]]$b, copula, tau)
  sh3 <- cif_copula_single(crn[[3]]$a, crn[[3]]$b, copula, tau)
  sh4 <- cif_copula_single(crn[[4]]$a, crn[[4]]$b, copula, tau)
  
  sh_l <- pmin(pmin(sh1, sh2), pmin(sh3, sh4))
  sh_u <- pmax(pmax(sh1, sh2), pmax(sh3, sh4))
  
  if (do_samp && B > 0) {
    
    for (b in seq_len(B)) {
      
      a <- runif(length(c1), c1_l, c1_u)
      d <- runif(length(c2), c2_l, c2_u)
      a <- trim_cif(a); d <- trim_cif(d)
      pr <- proj_sum(a, d)
      
      shb <- cif_copula_single(pr$a, pr$b, copula, tau)
      sh_l <- pmin(sh_l, shb)
      sh_u <- pmax(sh_u, shb)
    }
  }
  
  list(mean = sh_mu, lwr = sh_l, upr = sh_u)
}