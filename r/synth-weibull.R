
# =====================================================================
# Weibull proportional-hazards synthetic DGM. Source alongside
# synth-shared.r (and, optionally, synth-lognormal.r -- no clash,
# everything family-specific here is either uniquely named or an S3
# method dispatched on class(g) == "dgm_weibull", set by
# gen_surv_weibull() below.
#
# Parameterization. Each of T1, T2, C is Weibull-PH with
#   eta   = mu0 + muX*X + Z%*%muZ + W%*%muW        (mk_mu, unchanged)
#   H(t)  = (t/scale)^shape * exp(eta),   S(t) = exp(-H(t))
# so `par` lists carry `shape` and `scale` where the log-normal file
# carried `sT`. Everything else (mu0/muX/muZ/muW/winsL) is identical.
#
# SIGN CONVENTION FLIPS relative to the AFT log-normal file: eta is a
# log-hazard-ratio, so a LARGER muX now means SHORTER survival / HIGHER
# incidence. Flip the signs of muX/muZ/muW when porting parameters.
#
# Why this family: with a shape shared across the two causes, the
# conditional CIF is available in closed form,
#   CIF_k(t | v) = (theta_k / Theta) * (1 - exp(-t^shape * Theta)),
#   theta_j = scale_j^{-shape} exp(eta_j),  Theta = theta_1 + theta_2,
# which (i) removes all quadrature error from the oracle, and (ii) keeps
# CIFs in the interior of (0,1) for moderate ||eta||, unlike the
# saturating log-normal race. Unequal shapes are supported via a
# Gauss-Legendre fallback (exact to ~1e-12), but the closed form is used
# whenever shapes match, and oracle_nu_cr REQUIRES matching shapes.
# =====================================================================


# ------- core curve helpers (weibull-specific) -------

# cumulative hazard and survival for a Weibull-PH margin
wb_H <- function(t, eta, par) outer(exp(eta), (t / par$scale)^par$shape, `*`)
wb_S <- function(t, eta, par) exp(-wb_H(t, eta, par))

# theta = scale^{-shape} * exp(eta): the "rate" multiplying t^shape
wb_theta <- function(eta, par) exp(eta) / par$scale^par$shape

weibull_survival_curves <- function(df, event_par, t_grid) {
  X <- if ("X" %in% names(df)) df$X else df$majority
  Z <- as.matrix(df[, grepl("^z\\d+$", names(df)), drop = FALSE])
  W <- as.matrix(df[, grepl("^w\\d+$", names(df)), drop = FALSE])

  eta <- mk_mu(X, Z, W, event_par)

  out <- wb_S(t_grid, eta, event_par)

  dimnames(out) <- list(NULL, paste0("t_", t_grid))
  out
}

# Backward-compatible wrapper for the original event time
S_T_curves_weibull <- function(df, par, t_grid) {
  weibull_survival_curves(df, par$T, t_grid)
}

# Latent cause-specific survival curves in the CR setting
S_T1_curves_weibull <- function(df, par, t_grid) {
  weibull_survival_curves(df, par$T, t_grid)
}

S_T2_curves_weibull <- function(df, par, t_grid) {
  stopifnot(!is.null(par$T2))
  weibull_survival_curves(df, par$T2, t_grid)
}

# All-cause event-free survival, excluding censoring
S_all_curves_weibull <- function(df, par, t_grid) {
  stopifnot(!is.null(par$T2))
  S_T1_curves_weibull(df, par, t_grid) *
    S_T2_curves_weibull(df, par, t_grid)
}

# Censoring survival
G_C_curves_weibull <- function(df, par, t_grid) {
  weibull_survival_curves(df, par$C, t_grid)
}

ground_truth_G.dgm_weibull <- function(g, t_grid) {
  G_C_curves_weibull(g$data, g$par, t_grid)
}


# ------- generator -------

gen_surv_weibull <- function(n, k=3, pz=3, q=2,
                     SigU = diag(k),
                     A = matrix(rnorm(pz*k), pz, k),
                     beta = rnorm(k),
                     alpha = rep(0.5, q),
                     B = matrix(rnorm(pz*q), pz, q),
                     sZ = rep(0.7, pz),
                     sW = rep(0.7, q),
                     # T params (Weibull-PH): eta = mu0 + muX*X + Z%*%muZ + W%*%muW
                     T_par = list(mu0=0, muX=-0.4, muZ=rnorm(pz,0,0.15),
                                  muW=rnorm(q,-0.15,0.15), shape=1.5, scale=1),
                     # C params (indep Weibull-PH) with different coefs
                     C_par = list(mu0=-0.4, muX=0.15, muZ=rnorm(pz,0,0.1),
                                  muW=rnorm(q,-0.1,0.1), shape=1.2, scale=1),
                     T2_par = NULL,
                     pW = NULL,
                     seed = NULL) {
  if(!is.null(seed)) set.seed(seed)

  # U ~ N_k(0, SigU)
  U <- MASS::mvrnorm(n, mu=rep(0,k), Sigma=SigU)

  # Z = A U + eps_Z
  EZ <- matrix(rnorm(n*pz, 0, rep(sZ, each=n)), n, pz)
  Z <- U %*% t(A) + EZ
  colnames(Z) <- paste0("z", seq_len(pz))

  # X = 1{ logit^{-1}(beta^T U) > U(0,1) }
  eta_x <- drop(U %*% beta / 3)
  p <- 1/(1+exp(-eta_x))
  X <- rbinom(n, 1, p)

  # W = B^T Z + eps_W + M*(alpha*X)  (mixture switch M)
  if (is.null(pW)) pW <- runif(1, 0.6, 0.75)
  MW <- rbinom(n, 1, pW)

  EW <- matrix(rnorm(n*q, 0, rep(sW, each=n)), n, q)
  W <- (Z %*% B) + EW + outer(X * MW, alpha)
  colnames(W) <- paste0("w", seq_len(q))

  # generate counterfactuals for W (reuse same MW and EW)
  Wx0 <- (Z %*% B) + EW
  Wx1 <- (Z %*% B) + EW + outer(MW, alpha)
  colnames(Wx0) <- colnames(Wx1) <- colnames(W)

  # T = scale * (E * exp(-eta))^(1/shape),  E ~ Exp(1)
  rwb <- function(eta, par) {
    par$scale * (rexp(length(eta)) * exp(-eta))^(1 / par$shape)
  }

  etaT1 <- mk_mu(X, Z, W, T_par)
  T1 <- rwb(etaT1, T_par)

  etaC <- mk_mu(X, Z, W, C_par)
  C <- rwb(etaC, C_par)

  if (is.null(T2_par)) {

    M <- pmin(T1, C)
    delta <- as.integer(T1 <= C)

    df <- data.frame(
      majority = X, Z, W,
      T = T1, C = C,
      event_time = M,
      event = delta
    )

  } else {

    etaT2 <- mk_mu(X, Z, W, T2_par)
    T2 <- rwb(etaT2, T2_par)

    M <- pmin(T1, T2, C)

    delta <- ifelse(
      C < pmin(T1, T2), 0L,
      ifelse(T1 <= T2, 1L, 2L)
    )

    df <- data.frame(
      majority = X, Z, W,
      T1 = T1, T2 = T2, C = C,
      event_time = M,
      event = delta
    )
  }

  out <- list(data=df,
       par=list(T=T_par, T2=T2_par, C=C_par, A=A, beta=beta, alpha=alpha, B=B,
                SigU=SigU, sZ=sZ, sW=sW, pW=pW),
       cfW = list(Wx0=Wx0, Wx1=Wx1))
  class(out) <- c("dgm_weibull", class(out))
  out
}


# ------- counterfactual data construction -------

potential_survival_curves_from_gen_weibull <- function(g, t_grid, event_par) {
  df <- g$data
  cfW <- g$cfW

  list(
    x0_wx0 = weibull_survival_curves(
      make_df_with(df, 0L, cfW$Wx0), event_par, t_grid
    ),
    x0_wx1 = weibull_survival_curves(
      make_df_with(df, 0L, cfW$Wx1), event_par, t_grid
    ),
    x1_wx0 = weibull_survival_curves(
      make_df_with(df, 1L, cfW$Wx0), event_par, t_grid
    ),
    x1_wx1 = weibull_survival_curves(
      make_df_with(df, 1L, cfW$Wx1), event_par, t_grid
    )
  )
}

S_T_potential_curves_from_gen.dgm_weibull <- function(g, t_grid) {
  out <- potential_survival_curves_from_gen_weibull(
    g,
    t_grid,
    g$par$T
  )

  names(out) <- c(
    "S_x0_wx0",
    "S_x0_wx1",
    "S_x1_wx0",
    "S_x1_wx1"
  )

  out
}


# ------- latent counterfactual event times for CR truth -------

# rank-preserving exogenous noise: E = H(T) ~ Exp(1)
latent_error_from_gen.dgm_weibull <- function(g, event) {
  stopifnot(event %in% c(1L, 2L))

  df <- g$data
  par_event <- if (event == 1L) g$par$T else g$par$T2
  time_col <- if (event == 1L) "T1" else "T2"

  stopifnot(
    !is.null(par_event),
    time_col %in% names(df)
  )

  X <- df$majority
  Z <- as.matrix(df[, grepl("^z\\d+$", names(df)), drop = FALSE])
  W <- as.matrix(df[, grepl("^w\\d+$", names(df)), drop = FALSE])

  eta_obs <- mk_mu(X, Z, W, par_event)

  (df[[time_col]] / par_event$scale)^par_event$shape * exp(eta_obs)
}

counterfactual_latent_time.dgm_weibull <- function(g, event, xy, wx,
                                       latent_error = NULL) {
  stopifnot(
    event %in% c(1L, 2L),
    xy %in% c(0L, 1L),
    wx %in% c(0L, 1L)
  )

  if (is.null(latent_error)) {
    latent_error <- latent_error_from_gen(g, event)
  }

  par_event <- if (event == 1L) g$par$T else g$par$T2
  W_cf <- g$cfW[[paste0("Wx", wx)]]

  df_cf <- make_df_with(
    g$data,
    X_new = xy,
    W_new = W_cf
  )

  X_cf <- df_cf$X
  Z <- as.matrix(
    df_cf[, grepl("^z\\d+$", names(df_cf)), drop = FALSE]
  )
  W <- as.matrix(
    df_cf[, grepl("^w\\d+$", names(df_cf)), drop = FALSE]
  )

  eta_cf <- mk_mu(X_cf, Z, W, par_event)

  par_event$scale * (latent_error * exp(-eta_cf))^(1 / par_event$shape)
}


# ------- oracle conditional CIF -------

# CIF_k(t | eta_e, eta_o) for Weibull-PH margins.
# Equal shapes  -> closed form (theta_e/Theta)(1 - exp(-t^p Theta)).
# Unequal shapes -> substitution u = t v^{1/p_e} gives
#   CIF_k(t) = H_e(t) * int_0^1 exp(-sum_j H_j(t v^{p_j/p_e})) dv,
# integrated by Gauss-Legendre on [0,1] (smooth, no endpoint singularity).
.wb_cif <- function(t_grid, eta_e, eta_o, par_e, par_o, n_gl = 64) {

  th_e <- wb_theta(eta_e, par_e)
  th_o <- wb_theta(eta_o, par_o)

  if (isTRUE(all.equal(par_e$shape, par_o$shape))) {

    p <- par_e$shape
    Th <- th_e + th_o
    ratio <- th_e / Th
    out <- outer(Th, t_grid^p, `*`)
    out <- sweep(-expm1(-out), 1, ratio, `*`)
    return(pmin(pmax(out, 0), 1))
  }

  gl <- statmod::gauss.quad(n_gl, kind = "legendre")
  nodes <- (gl$nodes + 1) / 2
  wts <- gl$weights / 2

  pe <- par_e$shape; po <- par_o$shape
  out <- matrix(0, length(eta_e), length(t_grid))

  for (tt in seq_along(t_grid)) {
    t <- t_grid[tt]
    acc <- numeric(length(eta_e))
    for (qq in seq_len(n_gl)) {
      v <- nodes[qq]
      He <- th_e * t^pe * v                # H_e(t v^{1/p_e}) = H_e(t) * v
      Ho <- th_o * t^po * v^(po / pe)      # H_o(t v^{1/p_e})
      acc <- acc + wts[qq] * exp(-(He + Ho))
    }
    out[, tt] <- th_e * t^pe * acc
  }

  pmin(pmax(out, 0), 1)
}

# Oracle conditional CIF: F_k(t | X:=xy, Z, W) using each unit's OWN
# (factual) Z, W -- only X is intervened.
CIF_conditional_exact.dgm_weibull <- function(g, t_grid, event, n_gl = 64, ...) {
  stopifnot(!is.null(g$par$T2), event %in% c(1L, 2L))

  par_e <- if (event == 1L) g$par$T else g$par$T2
  par_o <- if (event == 1L) g$par$T2 else g$par$T

  df <- g$data
  Z <- as.matrix(df[, grepl("^z\\d+$", names(df)), drop = FALSE])
  W <- as.matrix(df[, grepl("^w\\d+$", names(df)), drop = FALSE])

  make_cif <- function(xy) {
    X <- rep(xy, nrow(Z))
    .wb_cif(t_grid, mk_mu(X, Z, W, par_e), mk_mu(X, Z, W, par_o),
            par_e, par_o, n_gl)
  }

  list(cifx0 = make_cif(0L), cifx1 = make_cif(1L))
}

# same object; kept under the old name (s_grid_n ignored -- no grid needed)
CIF_conditional_from_gen.dgm_weibull <- function(g, t_grid, event, s_grid_n = NULL, ...) {
  CIF_conditional_exact(g, t_grid, event)
}


# ------- oracle nested CIF regression nu(xy,xw,z) -------

# nu = E[ CIF_k(t | xy, Z, W) | X = xw, Z ]. Only the pair of linear
# predictors (eta_e, eta_o) matters, and given Z (and the M_W mixture
# component) that pair is bivariate normal, so the W-integral is a 2-d
# Gauss-Hermite quadrature. Requires a shared shape across causes, which
# is what makes the inner CIF closed-form.
.oracle_cif_integrated_W_weibull <- function(Z, meanW, xy, par_e, par_o, sW,
                                     t_grid, n_gh = 12) {
  stopifnot(is.null(par_e$winsL), is.null(par_o$winsL), all(t_grid > 0),
            isTRUE(all.equal(par_e$shape, par_o$shape)))

  p <- par_e$shape

  # marginal means / covariance of (eta_e, eta_o) given Z, component mean
  me <- par_e$mu0 + par_e$muX * xy + drop(Z %*% par_e$muZ) +
    drop(meanW %*% par_e$muW)
  mo <- par_o$mu0 + par_o$muX * xy + drop(Z %*% par_o$muZ) +
    drop(meanW %*% par_o$muW)

  ve <- sum((par_e$muW * sW)^2)
  vo <- sum((par_o$muW * sW)^2)
  ceo <- sum(par_e$muW * par_o$muW * sW^2)

  se <- sqrt(ve)
  b_cond <- if (ve > 0) ceo / se else 0
  so_cond <- sqrt(max(vo - (if (ve > 0) ceo^2 / ve else 0), 0))

  gh <- statmod::gauss.quad(n_gh, kind = "hermite")
  nd <- sqrt(2) * gh$nodes
  wt <- gh$weights / sqrt(pi)

  tp <- t_grid^p
  sc_e <- par_e$scale^p
  sc_o <- par_o$scale^p

  out <- matrix(0, nrow(Z), length(t_grid))

  for (i in seq_len(n_gh)) {
    for (j in seq_len(n_gh)) {

      eta_e <- me + se * nd[i]
      eta_o <- mo + b_cond * nd[i] + so_cond * nd[j]

      th_e <- exp(eta_e) / sc_e
      th_o <- exp(eta_o) / sc_o
      Th <- th_e + th_o
      ratio <- th_e / Th

      cif <- sweep(-expm1(-outer(Th, tp, `*`)), 1, ratio, `*`)
      out <- out + (wt[i] * wt[j]) * cif
    }
  }

  pmin(pmax(out, 0), 1)
}

oracle_nu_cr.dgm_weibull <- function(g, t_grid, event, xw, xy, n_gl = 12, ...) {
  stopifnot(event %in% c(1L, 2L), xw %in% 0:1, xy %in% 0:1)

  par_e <- if (event == 1L) g$par$T else g$par$T2
  par_o <- if (event == 1L) g$par$T2 else g$par$T

  Z <- as.matrix(g$data[, grepl("^z\\d+$", names(g$data)), drop = FALSE])
  meanW0 <- Z %*% g$par$B
  meanW1 <- sweep(meanW0, 2, xw * g$par$alpha, `+`)

  # M_W = 0 / 1 mixture
  nu0 <- .oracle_cif_integrated_W_weibull(
    Z, meanW0, xy, par_e, par_o, g$par$sW, t_grid, n_gl
  )

  if (xw == 0L) return(nu0)

  nu1 <- .oracle_cif_integrated_W_weibull(
    Z, meanW1, xy, par_e, par_o, g$par$sW, t_grid, n_gl
  )

  (1 - g$par$pW) * nu0 + g$par$pW * nu1
}
