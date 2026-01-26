
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

chf_rfs_cf <- function(data, X, time_var, event_var, rhs, time_interest,
                       balance_groups, split_forest, K = 5, folds = NULL, ...) {
  n <- nrow(data)
  idx <- sample(rep(1:K, length.out = n))
  
  if (!is.null(folds)) {
    assert_that(length(folds) == K)
    val_indx1 <- vector("list", K)
  }
  
  if (balance_groups) {
    
    wt0 <- 1 / mean(data[[X]] == 0)
    wt1 <- 1 / mean(data[[X]] == 1) 
  }
  
  chf <- chfx0 <- chfx1 <- srv <- srvx0 <- srvx1 <- 
    cif <- cifx0 <- cifx1 <- vector("list", K)
  if (is.null(time_interest)) time_interest <- 150 # use 150 time points
  ind <- NULL
  
  # check if any competing risk events are considered
  is_cr <- if (!all(data[[event_var]] %in% c(0, 1))) TRUE else FALSE
  
  # instantiate the formula (rhs already aware of split_forest argument)
  frml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  
  for (k in seq_len(K)) {
    
    if (!is.null(folds)) {
      
      trn <- data[folds[[k]][["dev"]]]
      val <- data[c(which(folds[[k]][["val"]]), which(folds[[k]][["tst"]]))]
      val_indx1[[k]] <- val[[X]] == 1
    } else {
      
      trn <- data[idx != k]
      val <- data[idx == k]
      ind <- c(ind, which(idx == k))
    }
    indx1 <- trn[[X]] == 1
    
    valx0 <- copy(val)
    valx1 <- copy(val)
    valx0[[X]] <- 0
    valx1[[X]] <- 1
    
    if (split_forest) {
      
      objx1 <- rfsrc(frml, data = trn[indx1, ], ntime = time_interest,
                     samptype = "swr", ...)
      if (length(time_interest) == 1) time_interest <- objx1$time.interest
      objx0 <- rfsrc(frml, data = trn[!indx1, ], ntime = time_interest,
                     samptype = "swr", ...)
    } else {
      
      # training weights
      gwt <- if (balance_groups) ifelse(trn[[X]], wt1, wt0) else rep(1, nrow(trn))
      gsize <- if (balance_groups) 2 * min(table(trn[[X]])) else nrow(trn)
      objx0 <- objx1 <- rfsrc(frml, data = trn, ntime = time_interest, 
                              case.wt = gwt, samptype = "swr", sampsize = gsize,
                              ...)
      if (length(time_interest) == 1) time_interest <- objx1$time.interest
    }
    
    predsx0 <- predict(objx0, newdata = valx0)
    predsx1 <- predict(objx1, newdata = valx1)
    
    chfx0[[k]] <- predsx0$chf
    chfx1[[k]] <- predsx1$chf
    
    if (is_cr) {
      
      cifx0[[k]] <- predsx0$cif
      cifx1[[k]] <- predsx1$cif
      
      # for competing risks, overall survival S(t) = 1 - \sum CIF_j(t)
      srvx0[[k]] <- 1 - apply(cifx0[[k]], c(1, 2), sum)
      srvx1[[k]] <- 1 - apply(cifx1[[k]], c(1, 2), sum)
    } else {
      
      srvx0[[k]] <- predsx0$survival
      srvx1[[k]] <- predsx1$survival
    }
    
    grid_ordx1 <- match_grids(time_interest, c(-Inf, objx1$time.interest))
    grid_ordx0 <- match_grids(time_interest, c(-Inf, objx0$time.interest))
    
    if (is_cr) {
      
      zeros_dim <- dim(chfx1[[k]])
      zeros_dim[2] <- 1
      zeros <- array(0, dim = zeros_dim)
      
      chfx0[[k]] <- abind(zeros, chfx0[[k]], along = 2)[, grid_ordx0, ]
      chfx1[[k]] <- abind(zeros, chfx1[[k]], along = 2)[, grid_ordx1, ]
      
      srvx0[[k]] <- cbind(1, srvx0[[k]])[, grid_ordx0]
      srvx1[[k]] <- cbind(1, srvx1[[k]])[, grid_ordx1]
      
      cifx0[[k]] <- abind(zeros, cifx0[[k]], along = 2)[, grid_ordx0, ]
      cifx1[[k]] <- abind(zeros, cifx1[[k]], along = 2)[, grid_ordx1, ]
      
    } else {
      
      chfx0[[k]] <- cbind(0, chfx0[[k]])[, grid_ordx0]
      chfx1[[k]] <- cbind(0, chfx1[[k]])[, grid_ordx1]
      
      srvx0[[k]] <- cbind(1, srvx0[[k]])[, grid_ordx0]
      srvx1[[k]] <- cbind(1, srvx1[[k]])[, grid_ordx1]
    }
  }
  
  indx1 <- data[[X]] == 1
  
  
  if (!is.null(folds)) { # if I have folds, returning val/tst separately!
    
    srv <- cif <- ret <- list()
    
    for (k in seq_len(K)) {
      
      val_idx <- seq_len(sum(folds[[k]][["val"]]))
      
      srv[[k]] <- srvx0[[k]]
      srv[[k]][val_indx1[[k]], ] <- srvx1[[k]][val_indx1[[k]], ]
      
      ret[[k]] <- list(
        srv_val = srv[[k]][val_idx, ],
        srvx0_val = srvx0[[k]][val_idx, ],
        srvx1_val = srvx1[[k]][val_idx, ],
        srv_tst = srv[[k]][-val_idx, ],
        srvx0_tst = srvx0[[k]][-val_idx, ],
        srvx1_tst = srvx1[[k]][-val_idx, ],
        time_interest = time_interest
      )
      
      if (is_cr) {
        
        cif[[k]] <- cifx0[[k]]
        cif[[k]][val_indx1[[k]], ,] <- cifx1[[k]][val_indx1[[k]], ,]
        
        add_cr <- list(
          cif_val = cif[[k]][val_idx, ,],
          cifx0_val = cifx0[[k]][val_idx, ,],
          cifx1_val = cifx1[[k]][val_idx, ,],
          cif_tst = cif[[k]][-val_idx, ,],
          cifx0_tst = cifx0[[k]][-val_idx, ,],
          cifx1_tst = cifx1[[k]][-val_idx, ,]
        )
        
        ret[[k]] <- c(ret[[k]], add_cr)
      }
    }
  } else { # without folds, returning out-of-fold predictions of everything
    
    if (is_cr) {
      
      ret <- list(
        chf = NULL,
        chfx0 = do.call(abind, args = list(chfx0, along = 1))[order(ind), ,],
        chfx1 = do.call(abind, args = list(chfx1, along = 1))[order(ind), ,],
        srv = NULL,
        srvx0 = do.call(rbind, srvx0)[order(ind), ],
        srvx1 = do.call(rbind, srvx1)[order(ind), ],
        cif = NULL,
        cifx0 = do.call(abind, args = list(cifx0, along = 1))[order(ind), ,],
        cifx1 = do.call(abind, args = list(cifx1, along = 1))[order(ind), ,],
        time_interest = time_interest
      )
      
      ret$chf <- ret$chfx0
      ret$chf[indx1, ,] <- ret$chfx1[indx1, ,]
      
      ret$srv <- ret$srvx0
      ret$srv[indx1, ] <- ret$srvx1[indx1, ]
      
      ret$cif <- ret$cifx0
      ret$cif[indx1, ,] <- ret$cifx1[indx1, ,]
    } else {
      
      ret <- list(
        chf = NULL,
        chfx0 = do.call(rbind, chfx0)[order(ind), ],
        chfx1 = do.call(rbind, chfx1)[order(ind), ],
        srv = NULL,
        srvx0 = do.call(rbind, srvx0)[order(ind), ],
        srvx1 = do.call(rbind, srvx1)[order(ind), ],
        cif = NULL, cifx0 = NULL, cifx1 = NULL,
        time_interest = time_interest
      )
      
      ret$chf <- ret$chfx0
      ret$chf[indx1, ] <- ret$chfx1[indx1, ]
      
      ret$srv <- ret$srvx0
      ret$srv[indx1, ] <- ret$srvx1[indx1, ]
    } 
  }
  
  ret
}

