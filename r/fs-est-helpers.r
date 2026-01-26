
# --- helpers used by all estimators -----------------------------------------
.fs_build_Y <- function(dt, time_var, event_var, tgrid, J) {
  n <- nrow(dt); Tn <- length(tgrid)
  Ysurv <- outer(dt[[time_var]], tgrid, ">=") * 1L
  if (J > 1L) {
    Ycif <- array(0L, dim = c(n, Tn, J))
    for (j in seq_len(J)) {
      idx <- dt[[event_var]] == j
      if (any(idx)) {
        when <- dt[[time_var]][idx]
        Ycif[idx, , j] <- outer(when, tgrid, "<=") * 1L
      }
    }
  } else Ycif <- NULL
  list(Ysurv = Ysurv, Ycif = Ycif)
}

.fs_measures_dt <- function(rows) data.table::rbindlist(rows)

# PS building blocks (reuse your exact e_z/e_zw/p1)
px_zw_z <- function(e_z, e_zw, x) if (x == 0) (1 - e_zw)/(1 - e_z) else e_zw/e_z
px_z_base <- function(e_z, e_zw, p1, zx) {
  if (is.na(zx)) return(rep(1, length(e_z)))
  if (zx == 0) (1 - e_zw)/(1 - p1) else e_z/p1
}

#' * MODEL-BASED (uses pack$S only) *
fs_estimate_model <- function(pack) {
  
  dt <- pack$data; X <- pack$X
  tgrid <- pack$grid$t; J <- pack$grid$J
  e_z <- pack$ps$e_z; e_zw <- pack$ps$e_zw; p1 <- pack$ps$p1
  
  do_scale <- function(scale) {
    if (scale == "surv") {
      m0 <- pack$S$S_x0; m1 <- pack$S$S_x1; ns <- 1L
    } else {
      if (is.null(pack$S$CIF_x0)) return(NULL)
      m0 <- pack$S$CIF_x0; m1 <- pack$S$CIF_x1; ns <- J
    }
    rows <- list()
    for (ev in seq_len(ns)) {
      getcol <- function(M) if (ns == 1L) M else M[,,ev]
      fx0 <- getcol(m0); fx1 <- getcol(m1)
      
      # named weights (depend only on wx, zx)
      w_wx0     <- px_zw_z(e_z, e_zw, 0) * px_z_base(e_z, e_zw, p1, NA)
      w_wx1     <- px_zw_z(e_z, e_zw, 1) * px_z_base(e_z, e_zw, p1, NA)
      w_wx0_zx0 <- px_zw_z(e_z, e_zw, 0) * px_z_base(e_z, e_zw, p1, 0)
      w_wx1_zx0 <- px_zw_z(e_z, e_zw, 1) * px_z_base(e_z, e_zw, p1, 0)
      w_wx1_zx1 <- px_zw_z(e_z, e_zw, 1) * px_z_base(e_z, e_zw, p1, 1)
      
      eval_w <- function(w, M, i) sum(w * M[, i]) / sum(w)
      
      idx1 <- dt[[X]] == 1L
      f_x0 <- colMeans(fx0[!idx1, , drop = FALSE])
      f_x1 <- colMeans(fx1[ idx1, , drop = FALSE])
      
      for (i in seq_along(tgrid)) {
        
        # tv-level: marginal 
        tv    <- f_x1[i] - f_x0[i]
        
        # population-level
        nde   <- eval_w(w_wx0,     fx1, i) - eval_w(w_wx0,     fx0, i)
        nie   <- eval_w(w_wx0,     fx1, i) - eval_w(w_wx1,     fx1, i)
        nse   <- (eval_w(w_wx1_zx1, fx1, i) - eval_w(w_wx1,     fx1, i)) -
          (eval_w(w_wx0_zx0, fx0, i) - eval_w(w_wx0,     fx0, i))
        
        # x-specific level
        ctfde <- eval_w(w_wx0_zx0, fx1, i) - eval_w(w_wx0_zx0, fx0, i)
        ctfie <- eval_w(w_wx0_zx0, fx1, i) - eval_w(w_wx1_zx0, fx1, i)
        ctfse <- eval_w(w_wx1_zx0, fx1, i) - eval_w(w_wx1_zx1, fx1, i)
        
        rows <- c(rows, list(
          list(time=tgrid[i], scale=scale, event=ev, effect="tv",    value=tv),
          list(time=tgrid[i], scale=scale, event=ev, effect="nde",   value=nde),
          list(time=tgrid[i], scale=scale, event=ev, effect="nie",   value=nie),
          list(time=tgrid[i], scale=scale, event=ev, effect="nse",   value=nse),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfde", value=ctfde),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfie", value=ctfie),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfse", value=ctfse)
        ))
      }
    }
    data.table::rbindlist(rows)
  }
  
  out <- Filter(Negate(is.null), list(do_scale("surv"), do_scale("cif")))
  list(measures = data.table::rbindlist(out), grid = pack$grid, meta = list(method="model"))
}

#' * IPW (PS + IPCW; uses observed Y with weights) *
fs_estimate_ipw <- function(pack) {
  
  dt <- pack$data; X <- pack$X
  tgrid <- pack$grid$t; J <- pack$grid$J
  e_z <- pack$ps$e_z; e_zw <- pack$ps$e_zw
  G0 <- pmax(pack$G$G_x0, 1e-6)
  G1 <- pmax(pack$G$G_x1, 1e-6)
  
  do_scale <- function(scale) {
    if (scale == "surv") {
      Y <- outer(dt[[pack$time_var]], tgrid, ">=") * 1L
      ns <- 1L
    } else {
      if (is.null(pack$S$CIF_x0)) return(NULL)
      ns <- J
      Y <- array(0L, dim = c(nrow(dt), length(tgrid), J))
      for (j in seq_len(J)) {
        idx <- dt[[pack$event_var]] == j
        if (any(idx))
          Y[idx, , j] <- outer(dt[[pack$time_var]][idx], tgrid, "<=") * 1L
      }
    }
    rows <- list()
    
    for (ev in seq_len(ns)) {
      getcol <- function(M) if (ns == 1L) M else M[,,ev]
      Yev <- getcol(Y)
      
      # named weights (depend only on wx, zx)
      w_wx0     <- px_zw_z(e_z, e_zw, 0) * px_z_base(e_z, e_zw, pack$ps$p1, NA)
      w_wx1     <- px_zw_z(e_z, e_zw, 1) * px_z_base(e_z, e_zw, pack$ps$p1, NA)
      w_wx0_zx0 <- px_zw_z(e_z, e_zw, 0) * px_z_base(e_z, e_zw, pack$ps$p1, 0)
      w_wx1_zx0 <- px_zw_z(e_z, e_zw, 1) * px_z_base(e_z, e_zw, pack$ps$p1, 0)
      w_wx1_zx1 <- px_zw_z(e_z, e_zw, 1) * px_z_base(e_z, e_zw, pack$ps$p1, 1)
      
      w_yx0 <- 1 / (1 - e_zw)
      w_yx1 <- 1 / e_zw
      
      idx1 <- dt[[X]] == 1L
      
      eval_ipw <- function(x, w, G, i)
        sum((dt[[X]] == x) * w * (Yev[, i] / G[, i])) /
        sum((dt[[X]] == x) * w)
      
      for (i in seq_along(tgrid)) {
        
        # factual survival by group
        tv <- mean((Yev[, i] / G1[, i])[idx1]) - mean((Yev[, i] / G0[, i])[!idx1])
        
        # population-level
        nde   <- eval_ipw(1, w_yx1 * w_wx0, G1, i) - eval_ipw(0, w_yx0 * w_wx0, G0, i)
        nie   <- eval_ipw(1, w_yx1 * w_wx0, G1, i) - eval_ipw(1, w_yx1 * w_wx1, G1, i)
        nse   <- eval_ipw(1, w_yx1 * w_wx1_zx1, G1, i) - eval_ipw(1, w_yx1 * w_wx1, G1, i) +
          eval_ipw(0, w_yx0 * w_wx0, G0, i) - eval_ipw(0, w_yx0 * w_wx0_zx0, G0, i)
        
        # x-specific level
        ctfde <- eval_ipw(1, w_yx1 * w_wx0_zx0, G1, i) - eval_ipw(0, w_yx0 * w_wx0_zx0, G0, i)
        ctfie <- eval_ipw(1, w_yx1 * w_wx0_zx0, G1, i) - eval_ipw(1, w_yx1 * w_wx1_zx0, G1, i)
        ctfse <- eval_ipw(1, w_yx1 * w_wx1_zx0, G1, i) - eval_ipw(1, w_yx1 * w_wx1_zx1, G1, i)
        
        # if (tgrid[i] > 9) browser()
        
        rows <- c(rows, list(
          list(time=tgrid[i], scale=scale, event=ev, effect="tv",    value=tv),
          list(time=tgrid[i], scale=scale, event=ev, effect="nde",   value=nde),
          list(time=tgrid[i], scale=scale, event=ev, effect="nie",   value=nie),
          list(time=tgrid[i], scale=scale, event=ev, effect="nse",   value=nse),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfde", value=ctfde),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfie", value=ctfie),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfse", value=ctfse)
        ))
      }
    }
    data.table::rbindlist(rows)
  }
  
  out <- Filter(Negate(is.null), list(do_scale("surv"), do_scale("cif")))
  list(measures = data.table::rbindlist(out), grid = pack$grid, meta = list(method="ipw"))
}

#' * AIPW (DR for treatment + censoring) -> still to be implemented *
fs_estimate_aipw <- function(pack) {
  dt <- pack$data; X <- pack$X
  tgrid <- pack$grid$t; J <- pack$grid$J
  e_z <- pack$ps$e_z; e_zw <- pack$ps$e_zw; p1 <- pack$ps$p1
  
  Y <- .fs_build_Y(dt, pack$time_var, pack$event_var, tgrid, J)
  Gf <- ifelse(dt[[X]] == 1, pack$G$G_x1, pack$G$G_x0); Gf <- pmax(Gf, 1e-6)
  
  add_block <- function(scale) {
    if (scale == "surv") { Ymat <- Y$Ysurv; m0 <- pack$S$S_x0; m1 <- pack$S$S_x1; ns <- 1L }
    else {
      if (is.null(Y$Ycif)) return(NULL)
      Ymat <- Y$Ycif;      m0 <- pack$S$CIF_x0; m1 <- pack$S$CIF_x1; ns <- J
    }
    rows <- list()
    for (ev in seq_len(ns)) {
      getcol <- function(M) if (ns == 1L) M else M[,,ev]
      Yev <- getcol(Ymat); m0e <- getcol(m0); m1e <- getcol(m1)
      
      wzw <- function(wx) .fs_w_zw_over_z(e_z, e_zw, wx)
      azx <- function(zx) .fs_adj_zx(e_z, e_zw, p1, zx)
      
      est_po <- function(fx, wx, zx, i) {
        w_med <- wzw(wx) * azx(zx)
        sfx   <- .fs_sfx(dt[[X]], e_z, p1, fx)
        mfx   <- if (fx==0) m0e[, i] else m1e[, i]
        # stabilized AIPW (censoring-robust)
        term <- mfx + sfx * (Yev[, i] - mfx) / Gf[, i]
        sum(w_med * term) / sum(w_med)
      }
      
      f_x0 <- est_po(0, NA, NA, seq_along(tgrid))
      f_x1 <- est_po(1, NA, NA, seq_along(tgrid))
      
      for (i in seq_along(tgrid)) {
        tv    <- f_x1[i] - f_x0[i]
        nde   <- est_po(1,0,NA,i) - est_po(0,0,NA,i)
        nie   <- est_po(1,0,NA,i) - est_po(1,1,NA,i)
        nse   <- (est_po(1,1,1,i) - est_po(1,1,NA,i)) - (est_po(0,0,0,i) - est_po(0,0,NA,i))
        ctfde <- est_po(1,0,0,i) - est_po(0,0,0,i)
        ctfie <- est_po(1,0,0,i) - est_po(1,1,0,i)
        ctfse <- est_po(1,1,0,i) - est_po(1,1,1,i)
        
        rows <- c(rows, list(
          list(time=tgrid[i], scale=scale, event=ev, effect="tv",    value=tv),
          list(time=tgrid[i], scale=scale, event=ev, effect="nde",   value=nde),
          list(time=tgrid[i], scale=scale, event=ev, effect="nie",   value=nie),
          list(time=tgrid[i], scale=scale, event=ev, effect="nse",   value=nse),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfde", value=ctfde),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfie", value=ctfie),
          list(time=tgrid[i], scale=scale, event=ev, effect="ctfse", value=ctfse)
        ))
      }
    }
    .fs_measures_dt(rows)
  }
  
  out <- Filter(Negate(is.null), list(add_block("surv"), add_block("cif")))
  list(measures = data.table::rbindlist(out), grid = pack$grid, meta = list(method="aipw"))
}
