
# ------- core helpers -------
lnorm_S <- function(t, mu, s) 1 - pnorm((log(t) - mu)/s)

mk_mu <- function(X, Z, W, par){
  mu <- drop(par$mu0 + par$muX*X + Z %*% par$muZ + W %*% par$muW)
  if (!is.null(par[["winsL"]])) {
    mu <- pmax(mu, par$winsL)
  }
  mu
}

# ------- generator -------
gen_surv <- function(n, k=3, pz=3, q=2,
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
  
  list(data=df,
       par=list(T=T_par, T2=T2_par,C=C_par, A=A, beta=beta, alpha=alpha, B=B,
                SigU=SigU, sZ=sZ, sW=sW, pW=pW),
       cfW = list(Wx0=Wx0, Wx1=Wx1))
}

# ------- generic latent-time curves -------

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
S_T_curves <- function(df, par, t_grid) {
  lognormal_survival_curves(df, par$T, t_grid)
}

# Latent cause-specific survival curves in the CR setting
S_T1_curves <- function(df, par, t_grid) {
  lognormal_survival_curves(df, par$T, t_grid)
}

S_T2_curves <- function(df, par, t_grid) {
  stopifnot(!is.null(par$T2))
  lognormal_survival_curves(df, par$T2, t_grid)
}

# All-cause event-free survival, excluding censoring
S_all_curves <- function(df, par, t_grid) {
  stopifnot(!is.null(par$T2))
  S_T1_curves(df, par, t_grid) *
    S_T2_curves(df, par, t_grid)
}

# Censoring survival
G_C_curves <- function(df, par, t_grid) {
  lognormal_survival_curves(df, par$C, t_grid)
}

ground_truth_G <- function(g, t_grid) {
  G_C_curves(g$data, g$par, t_grid)
}


# ------- counterfactual data construction -------

make_df_with <- function(df, X_new, W_new) {
  df2 <- df
  df2$X <- as.integer(X_new)
  
  wcols <- grepl("^w\\d+$", names(df2))
  df2[, wcols] <- W_new
  df2
}

potential_survival_curves_from_gen <- function(g, t_grid, event_par) {
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

# Backward-compatible wrapper
S_T_potential_curves_from_gen <- function(g, t_grid) {
  out <- potential_survival_curves_from_gen(
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


# ------- generic effect calculation -------

ground_truth_effects <- function(y00, y01, y10, y11, idx1,
                                 t_grid, scale, event) {
  
  tv <- colMeans(y11[idx1, , drop = FALSE]) -
    colMeans(y00[!idx1, , drop = FALSE])
  
  nde <- colMeans(y10 - y00)
  nie <- colMeans(y10 - y11)
  
  nse <-
    colMeans(y11[idx1, , drop = FALSE]) -
    colMeans(y11) -
    (
      colMeans(y00[!idx1, , drop = FALSE]) -
        colMeans(y00)
    )
  
  ctfde <- colMeans(
    (y10 - y00)[!idx1, , drop = FALSE]
  )
  
  ctfie <- colMeans(
    (y10 - y11)[!idx1, , drop = FALSE]
  )
  
  ctfse <-
    colMeans(y11[!idx1, , drop = FALSE]) -
    colMeans(y11[idx1, , drop = FALSE])
  
  out <- rbind(
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "tv", value = tv
    ),
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "nde", value = nde
    ),
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "nie", value = nie
    ),
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "nse", value = nse
    ),
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "ctfde", value = ctfde
    ),
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "ctfie", value = ctfie
    ),
    data.frame(
      time_interest = t_grid, scale = scale, event = event,
      effect = "ctfse", value = ctfse
    )
  )
  
  rownames(out) <- NULL
  data.table::as.data.table(out)
}


# ------- NIC ground truth -------

ground_truth <- function(g, tg) {
  Scf <- S_T_potential_curves_from_gen(g, tg)
  idx1 <- g$data$majority == 1
  
  ground_truth_effects(
    y00 = Scf$S_x0_wx0,
    y01 = Scf$S_x0_wx1,
    y10 = Scf$S_x1_wx0,
    y11 = Scf$S_x1_wx1,
    idx1 = idx1,
    t_grid = tg,
    scale = "surv",
    event = 1L
  )
}

# raw potential-outcome cells S(t|xy, W_xw), pre-differencing, for the
# same purpose as ground_truth_cr_cells below but on the survival scale
ground_truth_cells <- function(g, tg) {
  Scf <- S_T_potential_curves_from_gen(g, tg)
  idx <- g$data$majority

  cells <- list()
  for (xz in c(0L, 1L)) for (xw in c(0L, 1L)) for (xy in c(0L, 1L)) {
    y <- Scf[[paste0("S_x", xy, "_wx", xw)]]
    cells[[length(cells) + 1]] <- data.table::data.table(
      xz = xz, xw = xw, xy = xy, time_interest = tg,
      truth = colMeans(y[idx == xz, , drop = FALSE])
    )
  }
  data.table::rbindlist(cells)
}


# ------- latent counterfactual event times for CR truth -------

latent_error_from_gen <- function(g, event) {
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

counterfactual_latent_time <- function(g, event, xy, wx,
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
CIF_conditional_from_gen <- function(g, t_grid, event, s_grid_n = 2000) {
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
CIF_conditional_exact <- function(g, t_grid, event, n_gl = 200) {
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

cif_matrix_from_latent_times <- function(T_event, T_other, t_grid) {
  wins_competition <- T_event <= T_other
  
  outer(
    T_event,
    t_grid,
    function(tt, t) as.numeric(tt <= t)
  ) * wins_competition
}

CIF_potential_curves_from_gen <- function(g, t_grid, event) {
  stopifnot(
    !is.null(g$par$T2),
    event %in% c(1L, 2L)
  )
  
  other_event <- if (event == 1L) 2L else 1L
  
  eps_event <- latent_error_from_gen(g, event)
  eps_other <- latent_error_from_gen(g, other_event)
  
  make_cif <- function(xy, wx) {
    T_event <- counterfactual_latent_time(
      g, event, xy, wx, eps_event
    )
    
    T_other <- counterfactual_latent_time(
      g, other_event, xy, wx, eps_other
    )
    
    cif_matrix_from_latent_times(
      T_event,
      T_other,
      t_grid
    )
  }
  
  list(
    CIF_x0_wx0 = make_cif(0L, 0L),
    CIF_x0_wx1 = make_cif(0L, 1L),
    CIF_x1_wx0 = make_cif(1L, 0L),
    CIF_x1_wx1 = make_cif(1L, 1L)
  )
}


# ------- competing-risks ground truth -------

ground_truth_cr <- function(g, tg, events = c(1L, 2L)) {
  stopifnot(!is.null(g$par$T2))
  
  idx1 <- g$data$majority == 1
  
  out <- lapply(events, function(event) {
    Ccf <- CIF_potential_curves_from_gen(
      g,
      tg,
      event = event
    )
    
    ground_truth_effects(
      y00 = Ccf$CIF_x0_wx0,
      y01 = Ccf$CIF_x0_wx1,
      y10 = Ccf$CIF_x1_wx0,
      y11 = Ccf$CIF_x1_wx1,
      idx1 = idx1,
      t_grid = tg,
      scale = "cif",
      event = event
    )
  })
  
  data.table::rbindlist(out)
}

# raw potential-outcome cells CIF_k(t|xy, W_xw), pre-differencing --
# same purpose as ground_truth_cr but skips the effect-differencing, so
# a specific psi(xz,xw,xy) cell can be checked in isolation
ground_truth_cr_cells <- function(g, tg, events = c(1L, 2L)) {
  stopifnot(!is.null(g$par$T2))

  idx <- g$data$majority

  out <- lapply(events, function(event) {
    Ccf <- CIF_potential_curves_from_gen(g, tg, event = event)

    cells <- list()
    for (xz in c(0L, 1L)) for (xw in c(0L, 1L)) for (xy in c(0L, 1L)) {
      y <- Ccf[[paste0("CIF_x", xy, "_wx", xw)]]
      cells[[length(cells) + 1]] <- data.table::data.table(
        xz = xz, xw = xw, xy = xy, event = event, time_interest = tg,
        truth = colMeans(y[idx == xz, , drop = FALSE])
      )
    }
    data.table::rbindlist(cells)
  })

  data.table::rbindlist(out)
}

# ------- oracle propensity scores -------

oracle_px_z <- function(g, n_gh = 40) {
  par <- g$par
  Z <- as.matrix(g$data[, grepl("^z\\d+$", names(g$data)), drop = FALSE])
  
  # U | Z is Gaussian
  SU <- par$SigU
  SZ <- par$A %*% SU %*% t(par$A) + diag(par$sZ^2)
  K <- SU %*% t(par$A) %*% solve(SZ)
  mu_uz <- Z %*% t(K)
  SU_z <- SU - K %*% par$A %*% SU
  SU_z <- (SU_z + t(SU_z)) / 2
  
  # eta = beta' U / 3 | Z is univariate Gaussian
  b <- par$beta / 3
  mu_eta <- drop(mu_uz %*% b)
  sd_eta <- sqrt(max(drop(t(b) %*% SU_z %*% b), 0))
  
  # E[expit(eta)] by Gauss-Hermite quadrature
  gh <- statmod::gauss.quad(n_gh, kind = "hermite")
  eta <- outer(mu_eta, sqrt(2) * sd_eta * gh$nodes, `+`)
  p1 <- drop(plogis(eta) %*% gh$weights) / sqrt(pi)
  
  list(1 - p1, p1)
}


.logspace_add <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

.row_logdnorm_diag <- function(x, mu, s) {
  z <- sweep(x - mu, 2, s, `/`)
  -0.5 * ncol(x) * log(2*pi) - sum(log(s)) - 0.5 * rowSums(z^2)
}

oracle_px_zw <- function(g, px_z = oracle_px_z(g)) {
  par <- g$par
  Z <- as.matrix(g$data[, grepl("^z\\d+$", names(g$data)), drop = FALSE])
  W <- as.matrix(g$data[, grepl("^w\\d+$", names(g$data)), drop = FALSE])
  
  mu0 <- Z %*% par$B
  
  # W | X=0,Z: mixture collapses to one Gaussian
  lf0 <- .row_logdnorm_diag(W, mu0, par$sW)
  
  # W | X=1,Z: (1-pW)N(mu0,SigmaW) + pW N(mu0+alpha,SigmaW)
  mu1 <- sweep(mu0, 2, par$alpha, `+`)
  lf1_shift <- .row_logdnorm_diag(W, mu1, par$sW)
  lf1 <- .logspace_add(log1p(-par$pW) + lf0,
                       log(par$pW) + lf1_shift)
  
  l0 <- log(px_z[[1]]) + lf0
  l1 <- log(px_z[[2]]) + lf1
  m <- pmax(l0, l1)
  p1 <- exp(l1 - m) / (exp(l0 - m) + exp(l1 - m))
  
  list(1 - p1, p1)
}


# ------- oracle nested CIF regression nu(xy,xw,z) -------

.oracle_cif_integrated_W <- function(Z, meanW, xy, par_e, par_o, sW,
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

oracle_nu_cr <- function(g, t_grid, event, xw, xy, n_gl = 80) {
  stopifnot(event %in% c(1L, 2L), xw %in% 0:1, xy %in% 0:1)
  
  par_e <- if (event == 1L) g$par$T else g$par$T2
  par_o <- if (event == 1L) g$par$T2 else g$par$T
  
  Z <- as.matrix(g$data[, grepl("^z\\d+$", names(g$data)), drop = FALSE])
  meanW0 <- Z %*% g$par$B
  meanW1 <- sweep(meanW0, 2, xw * g$par$alpha, `+`)
  
  # M_W = 0 / 1 mixture
  nu0 <- .oracle_cif_integrated_W(
    Z, meanW0, xy, par_e, par_o, g$par$sW, t_grid, n_gl
  )
  nu1 <- .oracle_cif_integrated_W(
    Z, meanW1, xy, par_e, par_o, g$par$sW, t_grid, n_gl
  )
  
  (1 - g$par$pW) * nu0 + g$par$pW * nu1
}


# convenience wrapper: same nesting as px_z, px_zw, ey_nest in cross_fit_surv
oracle_nuisances_cr <- function(g, t_grid, n_gh = 40, n_gl = 80) {
  px_z <- oracle_px_z(g, n_gh)
  px_zw <- oracle_px_zw(g, px_z)
  
  ey_nest <- lapply(0:1, function(xw)
    lapply(0:1, function(xy) {
      mats <- lapply(1:2, function(j)
        oracle_nu_cr(g, t_grid, event = j, xw = xw, xy = xy, n_gl = n_gl))
      lapply(seq_along(t_grid), function(t)
        lapply(1:2, function(j) mats[[j]][, t]))
    })
  )
  
  list(px_z = px_z, px_zw = px_zw, ey_nest = ey_nest)
}