
# ---------- xgboost implementation ----------
mm_from_rhs <- function(df, rhs) {
  # design matrix from rhs (no Surv() on LHS here)
  X <- model.matrix(as.formula(paste0("~", rhs)), data = df)
  # drop intercept if present
  if (colnames(X)[1] == "(Intercept)") X <- X[, -1, drop = FALSE]
  X
}

breslow_from_offset <- function(y, d, lp) {
  fit <- coxph(Surv(y, d) ~ offset(lp), ties = "breslow")
  basehaz(fit, centered = FALSE) # data.frame: time, hazard
}

step_eval <- function(x, x0, y0) {
  # right-continuous step function evaluation on grid x using (x0, y0)
  approx(x0, y0, xout = x, method = "constant", f = 0, rule = 2)$y
}

match_grids_xgb <- function(target, source) {
  # for each target time, index of the last source time <= target
  # source is assumed sorted; may contain -Inf as first element
  findInterval(target, source, left.open = FALSE, rightmost.closed = TRUE)
}

surv_chf_from_lp <- function(lp, H0_vec) {
  # returns list(chf=|T|x|n|, srv=|T|x|n|)
  # lp: length n, H0_vec: length T
  lam <- exp(lp)
  chf <- outer(lam, H0_vec)            # |T| x n
  srv <- exp(-chf)
  browser()
  list(chf = chf, srv = srv)
}

# cause-specific CIF via Aalen–Johansen (using cause-specific H0_k and lp_k)
cif_from_cs <- function(H0k_list, lp_list, tt) {
  K <- length(H0k_list)
  n <- length(lp_list[[1]])
  Lam_k <- lapply(1:K, function(k) outer(H0k_list[[k]], exp(lp_list[[k]])))  # |T| x n
  Lam_sum <- Reduce(`+`, Lam_k)
  S <- exp(-Lam_sum)                                    # |T| x n
  dH0 <- lapply(H0k_list, function(h) c(h[1], diff(h))) # |T|
  dLam_k <- lapply(1:K, function(k) outer(dH0[[k]], exp(lp_list[[k]])))
  S_lag <- rbind(rep(1, n), S[-nrow(S), , drop = FALSE])
  CIF_k <- lapply(1:K, function(k) apply(S_lag * dLam_k[[k]], 2, cumsum))    # each |T| x n
  list(CIF = CIF_k, S = S, tt = tt)
}

# build labels for XGBoost Cox with "negative time = censored"
xgb_cox_label <- function(time, event01) ifelse(event01 == 1, time, -time)

# fit one Cox-XGB model and return booster and linear predictors
fit_xgb_cox <- function(X_tr, y_tr, d_tr, params, nrounds, wt = NULL) {
  # early stopping controls (optional; live inside params)
  val_frac <- if (!is.null(params$val_frac)) params$val_frac else NULL
  esr      <- if (!is.null(params$early_stopping_rounds)) params$early_stopping_rounds else NULL
  seed     <- if (!is.null(params$seed)) params$seed else 1L
  params$val_frac <- params$early_stopping_rounds <- params$seed <- NULL  # not xgb params
  
  lab <- xgb_cox_label(y_tr, d_tr)
  if (is.null(wt)) wt <- rep(1, length(lab))
  
  # split only if early stopping requested
  if (!is.null(esr) && !is.null(val_frac) && val_frac > 0 && val_frac < 1) {
    set.seed(seed)
    n  <- nrow(X_tr)
    nv <- max(1L, floor(val_frac * n))
    idv <- sample.int(n, nv)
    
    dtr <- xgb.DMatrix(data = X_tr[-idv, , drop = FALSE],
                       label = lab[-idv], weight = wt[-idv])
    dva <- xgb.DMatrix(data = X_tr[ idv, , drop = FALSE],
                       label = lab[ idv], weight = wt[ idv])
    
    bst <- xgb.train(
      params = params,
      data   = dtr,
      nrounds = nrounds,
      watchlist = list(train = dtr, val = dva),
      early_stopping_rounds = esr,
      verbose = 1
    )
  } else {
    dtr <- xgb.DMatrix(data = X_tr, label = lab, weight = wt)
    bst <- xgb.train(params, dtr, nrounds = nrounds, verbose = 0)
  }
  
  list(bst = bst)
}

# predict linear predictor
pred_lp <- function(bst, X_new) as.numeric(predict(bst, X_new, outputmargin=TRUE))

# expand / align (CH, S, CIF) to a common time grid with initial (0 / 1) column
align_single <- function(mat, tt_src, tt_targ, is_chf = TRUE) {
  # mat: |T_src| x n  (time-major); returns n x |T_targ|
  # prepend 0 (CH) or 1 (S) at time -Inf like your RF code
  init <- if (is_chf) 0 else 1
  mat2 <- rbind(rep(init, ncol(mat)), mat)        # add -Inf row
  ord <- match_grids_xgb(tt_targ, c(-Inf, tt_src))
  t(mat2[ord, , drop = FALSE])                    # n x |T_targ|
}

align_cif <- function(arr, tt_src, tt_targ) {
  # arr: |T_src| x n x K  -> return n x |T_targ| x K ; prepend ones at t=-Inf per your code
  K <- dim(arr)[3]
  n <- dim(arr)[2]
  ones <- array(1, dim = c(1, n, K))
  arr2 <- abind(ones, arr, along = 1)
  ord <- match_grids_xgb(tt_targ, c(-Inf, tt_src))
  # subset time and permute to n x T x K
  a <- arr2[ord, , , drop = FALSE]
  a <- aperm(a, c(2, 1, 3))
  a
}

xgb_surv_cf <- function(data, X, time_var, event_var, rhs, time_interest,
                        balance_groups, split_forest, K = 5,
                        nrounds = 300,
                        params = list(objective = "survival:cox",
                                      eval_metric = "cox-nloglik",
                                      eta = 0.05, max_depth = 3,
                                      subsample = 0.8, colsample_bytree = 0.8,
                                      early_stopping_rounds=15, val_frac = 0.2),
                        ...) {
  n <- nrow(data)
  idx <- sample(rep(1:K, length.out = n))
  if (balance_groups) {
    wt0 <- 1 / mean(data[[X]] == 0)
    wt1 <- 1 / mean(data[[X]] == 1)
  }
  
  chf <- chfx0 <- chfx1 <- srv <- srvx0 <- srvx1 <-
    cif <- cifx0 <- cifx1 <- vector("list", K)
  if (is.null(time_interest)) time_interest <- 150
  ind <- NULL
  
  # competing risks?
  is_cr <- if (!all(data[[event_var]] %in% c(0, 1))) TRUE else FALSE
  
  # precompute design matrix builder over arbitrary data.frames
  # rhs is expected to include X term already
  for (k in 1:K) {
    trn <- data[idx != k, , drop = FALSE]
    val <- data[idx == k, , drop = FALSE]
    ind <- c(ind, which(idx == k))
    indx1 <- trn[[X]] == 1
    
    # counterfactual copies for validation
    valx0 <- copy(val); valx0[[X]] <- 0
    valx1 <- copy(val); valx1[[X]] <- 1
    
    # build matrices
    X_tr <- mm_from_rhs(trn, rhs)
    X_valx0 <- mm_from_rhs(valx0, rhs)
    X_valx1 <- mm_from_rhs(valx1, rhs)
    
    y_tr <- trn[[time_var]]
    d_tr_all <- trn[[event_var]]
    y_val <- val[[time_var]]
    
    # weights for pooled fit
    gwt <- if (balance_groups) ifelse(trn[[X]] == 1, wt1, wt0) else rep(1, nrow(trn))
    
    if (!is_cr) {
      # ------- single-risk case -------
      d_tr <- as.integer(d_tr_all == 1)
      
      if (split_forest) {
        # fit separate by treatment
        X_tr1 <- X_tr[indx1, , drop = FALSE]
        y_tr1 <- y_tr[indx1]
        d_tr1 <- d_tr[indx1]
        X_tr0 <- X_tr[!indx1, , drop = FALSE]
        y_tr0 <- y_tr[!indx1]
        d_tr0 <- d_tr[!indx1]
        
        fit1 <- fit_xgb_cox(X_tr1, y_tr1, d_tr1, params, nrounds)
        fit0 <- fit_xgb_cox(X_tr0, y_tr0, d_tr0, params, nrounds)
        
        lp_tr1 <- pred_lp(fit1$bst, X_tr1)
        lp_tr0 <- pred_lp(fit0$bst, X_tr0)
        
        # Breslow baselines
        bh1 <- breslow_from_offset(y_tr1, d_tr1, lp_tr1)
        bh0 <- breslow_from_offset(y_tr0, d_tr0, lp_tr0)
        
        # prediction lps for val counterfactuals
        lp_valx1 <- pred_lp(fit1$bst, X_valx1)
        lp_valx0 <- pred_lp(fit0$bst, X_valx0)
        
        # choose / create a time grid
        if (length(time_interest) == 1L) {
          
          #' * note: time interest could be handled better? *
          tt1 <- bh1$time; tt0 <- bh0$time
          tt <- seq(0, max(c(tt0, tt1)), length.out = time_interest)
        } else {
          tt <- time_interest
        }
        H01 <- step_eval(tt, bh1$time, bh1$hazard)
        H00 <- step_eval(tt, bh0$time, bh0$hazard)
        
        sc1 <- surv_chf_from_lp(lp_valx1, H01)
        sc0 <- surv_chf_from_lp(lp_valx0, H00)
        
        # align to "RF-like" shape: add initial column (0 or 1) and reorder dims -> n x |T|
        chfx1[[k]] <- align_single(sc1$chf, tt, tt, is_chf = TRUE)
        chfx0[[k]] <- align_single(sc0$chf, tt, tt, is_chf = TRUE)
        srvx1[[k]] <- align_single(sc1$srv, tt, tt, is_chf = FALSE)
        srvx0[[k]] <- align_single(sc0$srv, tt, tt, is_chf = FALSE)
        
        time_interest <- tt
        grid_ordx1 <- seq_along(tt) # already matched
        grid_ordx0 <- seq_along(tt)
      } else {
        # pooled fit
        fitp <- fit_xgb_cox(X_tr, y_tr, d_tr, params, nrounds,
                            wt = gwt)
        lp_tr <- pred_lp(fitp$bst, X_tr)
        bh <- breslow_from_offset(y_tr, d_tr, lp_tr)
        
        # time grid
        if (length(time_interest) == 1L) {
          tt <- seq(0, max(bh$time), length.out = time_interest)
        } else {
          tt <- time_interest
        }
        H0 <- step_eval(tt, bh$time, bh$hazard)
        
        # browser()
        
        plot(bh$time, bh$hazard, pch=19)
        lines(tt, H0, col = "blue")
        
        # lps for val under X=0 and X=1
        lp_valx0 <- pred_lp(fitp$bst, X_valx0)
        lp_valx1 <- pred_lp(fitp$bst, X_valx1)
        
        sc0 <- surv_chf_from_lp(lp_valx0, H0)
        sc1 <- surv_chf_from_lp(lp_valx1, H0)
        
        chfx0[[k]] <- align_single(sc0$chf, tt, tt, is_chf = TRUE)
        chfx1[[k]] <- align_single(sc1$chf, tt, tt, is_chf = TRUE)
        srvx0[[k]] <- align_single(sc0$srv, tt, tt, is_chf = FALSE)
        srvx1[[k]] <- align_single(sc1$srv, tt, tt, is_chf = FALSE)
        
        time_interest <- tt
        grid_ordx1 <- seq_along(tt)
        grid_ordx0 <- seq_along(tt)
      }
      
      # prepend initial states to mimic RF outputs already handled in align_single()
      
    } else {
      # ------- competing risks -------
      # define causes: assume event_var in {0,1,...,Kc}, with 0=censor
      causes <- sort(setdiff(unique(d_tr_all), 0L))
      Kc <- length(causes)
      
      # per-cause training indicators
      d_k_list <- lapply(causes, function(kc) as.integer(d_tr_all == kc))
      
      if (split_forest) {
        # fit per treatment, per cause
        # X=1
        X_tr1 <- X_tr[indx1, , drop = FALSE]; y_tr1 <- y_tr[indx1]
        X_valx1 <- X_valx1
        fits1 <- lapply(d_k_list, function(dk) {
          dk1 <- dk[indx1]
          fit_xgb_cox(X_tr1, y_tr1, dk1, params, nrounds)
        })
        lp_tr1_list <- lapply(seq_along(fits1), function(j) pred_lp(fits1[[j]]$bst, X_tr1))
        bh1_list <- mapply(function(lp, dk) breslow_from_offset(y_tr1, dk[indx1], lp),
                           lp_tr1_list, d_k_list, SIMPLIFY = FALSE)
        lp_valx1_list <- lapply(seq_along(fits1), function(j) pred_lp(fits1[[j]]$bst, X_valx1))
        
        # X=0
        X_tr0 <- X_tr[!indx1, , drop = FALSE]; y_tr0 <- y_tr[!indx1]
        X_valx0 <- X_valx0
        fits0 <- lapply(d_k_list, function(dk) {
          dk0 <- dk[!indx1]
          fit_xgb_cox(X_tr0, y_tr0, dk0, params, nrounds)
        })
        lp_tr0_list <- lapply(seq_along(fits0), function(j) pred_lp(fits0[[j]]$bst, X_tr0))
        bh0_list <- mapply(function(lp, dk) breslow_from_offset(y_tr0, dk[!indx1], lp),
                           lp_tr0_list, d_k_list, SIMPLIFY = FALSE)
        lp_valx0_list <- lapply(seq_along(fits0), function(j) pred_lp(fits0[[j]]$bst, X_valx0))
        
        # common time grid
        if (length(time_interest) == 1L) {
          tt_all <- sort(unique(unlist(c(lapply(bh1_list, `[[`, "time"),
                                         lapply(bh0_list, `[[`, "time")))))
          tt <- seq(0, max(tt_all), length.out = time_interest)
        } else {
          tt <- time_interest
        }
        H01_list <- lapply(bh1_list, function(bh) step_eval(tt, bh$time, bh$hazard))
        H00_list <- lapply(bh0_list, function(bh) step_eval(tt, bh$time, bh$hazard))
        
        # CIFs + S
        res1 <- cif_from_cs(H01_list, lp_valx1_list, tt)
        res0 <- cif_from_cs(H00_list, lp_valx0_list, tt)
        
        # shapes to match your RF post-processing
        # chf per cause: cumulative hazards per cause (needed downstream)
        chf1_arr <- abind(lapply(1:Kc, function(j) t(outer(H01_list[[j]], exp(lp_valx1_list[[j]])))), along = 3)
        chf0_arr <- abind(lapply(1:Kc, function(j) t(outer(H00_list[[j]], exp(lp_valx0_list[[j]])))), along = 3)
        # align (prepend zero plane at t=-Inf)
        chfx1[[k]] <- align_cif(aperm(chf1_arr, c(2,1,3)), tt, tt)  # n x T x Kc
        chfx0[[k]] <- align_cif(aperm(chf0_arr, c(2,1,3)), tt, tt)
        
        # srv (overall)
        srvx1[[k]] <- align_single(res1$S, tt, tt, is_chf = FALSE)
        srvx0[[k]] <- align_single(res0$S, tt, tt, is_chf = FALSE)
        
        # cif arrays per cause
        # res$CIF[[j]] is |T| x n ; stack to n x |T| x Kc
        cif1_arr <- abind(lapply(res1$CIF, function(M) t(M)), along = 3)
        cif0_arr <- abind(lapply(res0$CIF, function(M) t(M)), along = 3)
        cifx1[[k]] <- align_cif(aperm(cif1_arr, c(2,1,3)), tt, tt)  # n x T x Kc (with leading 1s at -Inf)
        cifx0[[k]] <- align_cif(aperm(cif0_arr, c(2,1,3)), tt, tt)
        
        time_interest <- tt
        grid_ordx1 <- seq_along(tt); grid_ordx0 <- seq_along(tt)
      } else {
        # pooled per-cause fits
        fits <- lapply(d_k_list, function(dk) fit_xgb_cox(X_tr, y_tr, dk, params, nrounds, wt = gwt))
        lp_tr_list <- lapply(seq_along(fits), function(j) pred_lp(fits[[j]]$bst, X_tr))
        bh_list <- mapply(function(lp, dk) breslow_from_offset(y_tr, dk, lp),
                          lp_tr_list, d_k_list, SIMPLIFY = FALSE)
        
        if (length(time_interest) == 1L) {
          tt_all <- sort(unique(unlist(lapply(bh_list, `[[`, "time"))))
          tt <- seq(0, max(tt_all), length.out = time_interest)
        } else {
          tt <- time_interest
        }
        H0_list <- lapply(bh_list, function(bh) step_eval(tt, bh$time, bh$hazard))
        
        # lps for val under X=0 and X=1 (same pooled fits)
        lp_valx0_list <- lapply(seq_along(fits), function(j) pred_lp(fits[[j]]$bst, X_valx0))
        lp_valx1_list <- lapply(seq_along(fits), function(j) pred_lp(fits[[j]]$bst, X_valx1))
        
        res0 <- cif_from_cs(H0_list, lp_valx0_list, tt)
        res1 <- cif_from_cs(H0_list, lp_valx1_list, tt)
        
        # CHFs per cause (needed)
        chf0_arr <- abind(lapply(1:length(H0_list), function(j) t(outer(H0_list[[j]], exp(lp_valx0_list[[j]])))), along = 3)
        chf1_arr <- abind(lapply(1:length(H0_list), function(j) t(outer(H0_list[[j]], exp(lp_valx1_list[[j]])))), along = 3)
        
        chfx0[[k]] <- align_cif(aperm(chf0_arr, c(2,1,3)), tt, tt)
        chfx1[[k]] <- align_cif(aperm(chf1_arr, c(2,1,3)), tt, tt)
        
        srvx0[[k]] <- align_single(res0$S, tt, tt, is_chf = FALSE)
        srvx1[[k]] <- align_single(res1$S, tt, tt, is_chf = FALSE)
        
        cif0_arr <- abind(lapply(res0$CIF, function(M) t(M)), along = 3)
        cif1_arr <- abind(lapply(res1$CIF, function(M) t(M)), along = 3)
        cifx0[[k]] <- align_cif(aperm(cif0_arr, c(2,1,3)), tt, tt)
        cifx1[[k]] <- align_cif(aperm(cif1_arr, c(2,1,3)), tt, tt)
        
        time_interest <- tt
        grid_ordx1 <- seq_along(tt); grid_ordx0 <- seq_along(tt)
      }
    }
    
    # note: RF version pads via 'grid_ord'—we aligned directly to 'tt'
    # retained grid_ord* variables for parity with your code, though unused now
  } # end CV loop
  
  indx1 <- data[[X]] == 1
  if (is_cr) {
    ret <- list(
      chf = NULL,
      chfx0 = do.call(abind, args = list(chfx0, along = 1))[order(ind), , ],
      chfx1 = do.call(abind, args = list(chfx1, along = 1))[order(ind), , ],
      srv = NULL,
      srvx0 = do.call(rbind, srvx0)[order(ind), ],
      srvx1 = do.call(rbind, srvx1)[order(ind), ],
      cif = NULL,
      cifx0 = do.call(abind, args = list(cifx0, along = 1))[order(ind), , ],
      cifx1 = do.call(abind, args = list(cifx1, along = 1))[order(ind), , ],
      time_interest = time_interest
    )
    ret$chf <- ret$chfx0
    ret$chf[indx1, , ] <- ret$chfx1[indx1, , ]
    ret$srv <- ret$srvx0
    ret$srv[indx1, ] <- ret$srvx1[indx1, ]
    ret$cif <- ret$cifx0
    ret$cif[indx1, , ] <- ret$cifx1[indx1, , ]
  } else {
    ret <- list(
      chf = NULL,
      chfx0 = do.call(rbind, chfx0)[order(ind), , drop = FALSE],
      chfx1 = do.call(rbind, chfx1)[order(ind), , drop = FALSE],
      srv = NULL,
      srvx0 = do.call(rbind, srvx0)[order(ind), , drop = FALSE],
      srvx1 = do.call(rbind, srvx1)[order(ind), , drop = FALSE],
      cif = NULL, cifx0 = NULL, cifx1 = NULL,
      time_interest = time_interest
    )
    ret$chf <- ret$chfx0
    ret$chf[indx1, ] <- ret$chfx1[indx1, ]
    ret$srv <- ret$srvx0
    ret$srv[indx1, ] <- ret$srvx1[indx1, ]
  }
  ret
}
