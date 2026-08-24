
# =====================================================================
# Log-normal (AFT) synthetic DGM. Source alongside synth-shared.r (and,
# optionally, synth-weibull.R -- no clash, everything family-specific
# here is either uniquely named or an S3 method dispatched on
# class(g) == "dgm_lognormal", set by gen_surv_lognormal() below.
# =====================================================================


# ------- core curve helpers (lognormal-specific) -------

lnorm_S <- function(t, mu, s) 1 - pnorm((log(t) - mu)/s)

lognormal_survival_curves <- function(df, event_par, t_grid) {
  X <- if ("X" %in% names(df)) df$X else df$majority
  Z <- as.matrix(df[, grepl("^z\\d+$", names(df)), drop = FALSE])
  W <- as.matrix(df[, grepl("^w\\d+$", names(df)), drop = FALSE])

  mu <- mk_mu(X, Z, W, event_par)

  out <- outer(
    mu,
    t_grid,
    function(m, t) lnorm_S(t, m, event_par$sT)
  )

  dimnames(out) <- list(NULL, paste0("t_", t_grid))
  out
}

# Backward-compatible wrapper for the original event time
S_T_curves_lognormal <- function(df, par, t_grid) {
  lognormal_survival_curves(df, par$T, t_grid)
}

# Latent cause-specific survival curves in the CR setting
S_T1_curves_lognormal <- function(df, par, t_grid) {
  lognormal_survival_curves(df, par$T, t_grid)
}

S_T2_curves_lognormal <- function(df, par, t_grid) {
  stopifnot(!is.null(par$T2))
  lognormal_survival_curves(df, par$T2, t_grid)
}

# All-cause event-free survival, excluding censoring
S_all_curves_lognormal <- function(df, par, t_grid) {
  stopifnot(!is.null(par$T2))
  S_T1_curves_lognormal(df, par, t_grid) *
    S_T2_curves_lognormal(df, par, t_grid)
}

# Censoring survival
G_C_curves_lognormal <- function(df, par, t_grid) {
  lognormal_survival_curves(df, par$C, t_grid)
}

ground_truth_G.dgm_lognormal <- function(g, t_grid) {
  G_C_curves_lognormal(g$data, g$par, t_grid)
}


# ------- generator -------

gen_surv_lognormal <- function(n, k=3, pz=3, q=2,
                     SigU = diag(k),
                     A = matrix(rnorm(pz*k), pz, k),
                     beta = rnorm(k),
                     alpha = rep(0.5, q),
                     B = matrix(rnorm(pz*q), pz, q),
                     sZ = rep(0.7, pz),
                     sW = rep(0.7, q),
                     # T params (log-normal): mu = mu0 + muX*X + Z%*%muZ + muW*W ; sd = sT (const)
                     T_par = list(mu0=0, muX=0.6, muZ=rnorm(pz,0,0.2),
                                  muW=rnorm(q,0.2,0.2), sT=0.6),
                     # C params (indep log-normal) with different coefs
                     C_par = list(mu0=0.3, muX=-0.2, muZ=rnorm(pz,0,0.2),
                                  muW=rnorm(q,0.2,0.2), sT=1, winsL = 0),
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

  muT1 <- mk_mu(X, Z, W, T_par)
  T1 <- rlnorm(n, meanlog = muT1, sdlog = T_par$sT)

  muC <- mk_mu(X, Z, W, C_par)
  C <- rlnorm(n, meanlog = muC, sdlog = C_par$sT)

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

    muT2 <- mk_mu(X, Z, W, T2_par)
    T2 <- rlnorm(n, meanlog = muT2, sdlog = T2_par$sT)

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
       par=list(T=T_par, T2=T2_par,C=C_par, A=A, beta=beta, alpha=alpha, B=B,
                SigU=SigU, sZ=sZ, sW=sW, pW=pW),
       cfW = list(Wx0=Wx0, Wx1=Wx1))
  class(out) <- c("dgm_lognormal", class(out))
  out
}

# legacy alias: existing scripts call gen_surv() expecting log-normal
gen_surv <- gen_surv_lognormal


# ------- counterfactual data construction -------

potential_survival_curves_from_gen_lognormal <- function(g, t_grid, event_par) {
  df <- g$data
  cfW <- g$cfW

  list(
    x0_wx0 = lognormal_survival_curves(
      make_df_with(df, 0L, cfW$Wx0), event_par, t_grid
    ),
    x0_wx1 = lognormal_survival_curves(
      make_df_with(df, 0L, cfW$Wx1), event_par, t_grid
    ),
    x1_wx0 = lognormal_survival_curves(
      make_df_with(df, 1L, cfW$Wx0), event_par, t_grid
    ),
    x1_wx1 = lognormal_survival_curves(
      make_df_with(df, 1L, cfW$Wx1), event_par, t_grid
    )
  )
}

S_T_potential_curves_from_gen.dgm_lognormal <- function(g, t_grid) {
  out <- potential_survival_curves_from_gen_lognormal(
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

latent_error_from_gen.dgm_lognormal <- function(g, event) {
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

  mu_obs <- mk_mu(X, Z, W, par_event)

  (log(df[[time_col]]) - mu_obs) / par_event$sT
}

counterfactual_latent_time.dgm_lognormal <- function(g, event, xy, wx,
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

  mu_cf <- mk_mu(X_cf, Z, W, par_event)

  exp(mu_cf + par_event$sT * latent_error)
}

# Oracle conditional CIF: F_k(t | X:=xy, Z, W) using each unit's OWN
# (factual) Z, W -- only X is intervened. This is the true conditional
# probability (race between independent log-normal latent times), NOT a
# per-unit realization -- integrates out each unit's own exogenous noise,
# via F_k(t|X,Z,W) = int_0^t f_event(s|X,Z,W) * S_other(s|X,Z,W) ds.
# (counterfactual_latent_time above uses rank-preserved individual noise,
# which is correct for population-averaged ground truth but NOT for a
# per-unit conditional-probability nuisance.)
CIF_conditional_from_gen.dgm_lognormal <- function(g, t_grid, event, s_grid_n = 2000, ...) {
  stopifnot(!is.null(g$par$T2), event %in% c(1L, 2L))

  par_event <- if (event == 1L) g$par$T else g$par$T2
  par_other <- if (event == 1L) g$par$T2 else g$par$T

  df <- g$data
  Z <- as.matrix(df[, grepl("^z\\d+$", names(df)), drop = FALSE])
  W <- as.matrix(df[, grepl("^w\\d+$", names(df)), drop = FALSE])
  n <- nrow(df)

  s_grid <- seq(1e-6, max(t_grid), length.out = s_grid_n)
  d_s <- diff(s_grid)

  make_cif <- function(xy) {
    X <- rep(xy, n)
    mu_event <- mk_mu(X, Z, W, par_event)
    mu_other <- mk_mu(X, Z, W, par_other)

    dens <- outer(mu_event, s_grid, function(mu, s) dlnorm(s, mu, par_event$sT))
    surv_other <- outer(mu_other, s_grid,
                        function(mu, s) 1 - plnorm(s, mu, par_other$sT))
    integrand <- dens * surv_other

    avg <- (integrand[, -1, drop = FALSE] +
            integrand[, -ncol(integrand), drop = FALSE]) / 2
    step <- sweep(avg, 2, d_s, `*`)
    cum <- cbind(0, t(apply(step, 1, cumsum)))

    t(apply(cum, 1, function(row) approx(s_grid, row, xout = t_grid)$y))
  }

  list(cifx0 = make_cif(0L), cifx1 = make_cif(1L))
}

# exact oracle CIF via u = (log s - mu_e)/sT substitution + Gauss-Legendre
CIF_conditional_exact.dgm_lognormal <- function(g, t_grid, event, n_gl = 200, ...) {
  par_e <- if (event == 1L) g$par$T else g$par$T2
  par_o <- if (event == 1L) g$par$T2 else g$par$T
  Z <- as.matrix(g$data[, grepl("^z\\d+$", names(g$data)), drop = FALSE])
  W <- as.matrix(g$data[, grepl("^w\\d+$", names(g$data)), drop = FALSE])
  gl <- statmod::gauss.quad(n_gl, kind = "legendre")
  make_cif <- function(xy) {
    mu_e <- mk_mu(rep(xy, nrow(Z)), Z, W, par_e)
    mu_o <- mk_mu(rep(xy, nrow(Z)), Z, W, par_o)
    sapply(t_grid, function(t) {
      up <- (log(t) - mu_e) / par_e$sT
      lo <- pmin(up, -8); mid <- (up + lo) / 2; half <- (up - lo) / 2
      val <- 0
      for (q in seq_len(n_gl)) {
        u <- mid + half * gl$nodes[q]
        val <- val + gl$weights[q] * dnorm(u) *
          pnorm((mu_e + par_e$sT * u - mu_o) / par_o$sT, lower.tail = FALSE)
      }
      pmin(pmax(half * val, 0), 1)
    })
  }
  list(cifx0 = make_cif(0L), cifx1 = make_cif(1L))
}


# ------- oracle nested CIF regression nu(xy,xw,z) -------

.oracle_cif_integrated_W_lognormal <- function(Z, meanW, xy, par_e, par_o, sW,
                                     t_grid, n_gl = 80) {
  stopifnot(is.null(par_e$winsL), is.null(par_o$winsL), all(t_grid > 0))

  # Marginal log-times after integrating Gaussian W noise
  me <- par_e$mu0 + par_e$muX * xy + drop(Z %*% par_e$muZ) +
    drop(meanW %*% par_e$muW)
  mo <- par_o$mu0 + par_o$muX * xy + drop(Z %*% par_o$muZ) +
    drop(meanW %*% par_o$muW)

  ve <- par_e$sT^2 + sum((par_e$muW * sW)^2)
  vo <- par_o$sT^2 + sum((par_o$muW * sW)^2)
  ceo <- sum(par_e$muW * par_o$muW * sW^2)

  se <- sqrt(ve)
  beta_cond <- ceo / ve
  so_cond <- sqrt(vo - ceo^2 / ve)

  gl <- statmod::gauss.quad(n_gl, kind = "legendre")
  out <- matrix(0, nrow(Z), length(t_grid))

  # CIF = int phi(u) P(L_other >= L_event | L_event) du
  for (tt in seq_along(t_grid)) {
    up <- (log(t_grid[tt]) - me) / se
    hi <- pmax(-8, pmin(up, 8))
    lo <- rep(-8, length(up))
    mid <- (hi + lo) / 2
    half <- (hi - lo) / 2

    acc <- numeric(length(up))
    for (q in seq_len(n_gl)) {
      u <- mid + half * gl$nodes[q]
      le <- me + se * u
      mo_cond <- mo + beta_cond * (le - me)
      acc <- acc + gl$weights[q] * dnorm(u) *
        pnorm((mo_cond - le) / so_cond)
    }
    out[, tt] <- half * acc
  }

  pmin(pmax(out, 0), 1)
}

oracle_nu_cr.dgm_lognormal <- function(g, t_grid, event, xw, xy, n_gl = 80, ...) {
  stopifnot(event %in% c(1L, 2L), xw %in% 0:1, xy %in% 0:1)

  par_e <- if (event == 1L) g$par$T else g$par$T2
  par_o <- if (event == 1L) g$par$T2 else g$par$T

  Z <- as.matrix(g$data[, grepl("^z\\d+$", names(g$data)), drop = FALSE])
  meanW0 <- Z %*% g$par$B
  meanW1 <- sweep(meanW0, 2, xw * g$par$alpha, `+`)

  # M_W = 0 / 1 mixture
  nu0 <- .oracle_cif_integrated_W_lognormal(
    Z, meanW0, xy, par_e, par_o, g$par$sW, t_grid, n_gl
  )
  nu1 <- .oracle_cif_integrated_W_lognormal(
    Z, meanW1, xy, par_e, par_o, g$par$sW, t_grid, n_gl
  )

  (1 - g$par$pW) * nu0 + g$par$pW * nu1
}
