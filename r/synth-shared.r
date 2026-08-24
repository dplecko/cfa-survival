
# =====================================================================
# Family-agnostic scaffolding, shared by synth-lognormal.r and
# synth-weibull.R. Source both family files alongside this one -- they
# no longer clash, since every function that differs by family is an
# S3 generic dispatched on class(g) (set by gen_surv_lognormal() /
# gen_surv_weibull() as "dgm_lognormal" / "dgm_weibull").
#
# Generic stubs live here; each family file supplies the methods.
# =====================================================================

# ------- S3 generics (methods in synth-lognormal.r / synth-weibull.R) -------

S_T_potential_curves_from_gen <- function(g, t_grid) UseMethod("S_T_potential_curves_from_gen")
ground_truth_G <- function(g, t_grid) UseMethod("ground_truth_G")
CIF_conditional_exact <- function(g, t_grid, event, ...) UseMethod("CIF_conditional_exact")
CIF_conditional_from_gen <- function(g, t_grid, event, ...) UseMethod("CIF_conditional_from_gen")
latent_error_from_gen <- function(g, event) UseMethod("latent_error_from_gen")
counterfactual_latent_time <- function(g, event, xy, wx, latent_error = NULL) UseMethod("counterfactual_latent_time")
oracle_nu_cr <- function(g, t_grid, event, xw, xy, ...) UseMethod("oracle_nu_cr")


# ------- core helpers (family-independent) -------

mk_mu <- function(X, Z, W, par){
  mu <- drop(par$mu0 + par$muX*X + Z %*% par$muZ + W %*% par$muW)
  if (!is.null(par[["winsL"]])) {
    mu <- pmax(mu, par$winsL)
  }
  mu
}

make_df_with <- function(df, X_new, W_new) {
  df2 <- df
  df2$X <- as.integer(X_new)

  wcols <- grepl("^w\\d+$", names(df2))
  df2[, wcols] <- W_new
  df2
}

cif_matrix_from_latent_times <- function(T_event, T_other, t_grid) {
  wins_competition <- T_event <= T_other

  outer(
    T_event,
    t_grid,
    function(tt, t) as.numeric(tt <= t)
  ) * wins_competition
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


# ------- NIC ground truth (survival scale) -------

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


# ------- competing-risks ground truth (built on the dispatched leaves) -------

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


# ------- oracle propensity scores (family-independent: only depend on
#         the U/Z/X/W structure, never on T/C/T2) -------

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


# ------- DGP acceptance check (family-independent: dispatches inside) -------

# Reject parameterizations whose conditional CIFs saturate at 0/1: those
# are the ones the forest cannot learn, and where E2 is driven by a few
# extreme units. Aim for most mass interior at the evaluation times.
cif_interiority_check <- function(g, t_grid, events = c(1L, 2L),
                                  probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  out <- lapply(events, function(k) {
    cf <- CIF_conditional_exact(g, t_grid, event = k)
    do.call(rbind, lapply(seq_along(t_grid), function(tt) {
      v <- c(cf$cifx0[, tt], cf$cifx1[, tt])
      data.frame(event = k, time_interest = t_grid[tt],
                 frac_lt_01 = mean(v < 0.01), frac_gt_60 = mean(v > 0.60),
                 t(setNames(quantile(v, probs), paste0("q", probs * 100))))
    }))
  })
  do.call(rbind, out)
}
