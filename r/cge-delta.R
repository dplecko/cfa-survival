# ---- Archimedean generators ------------------------------------------------
# Returns list(phi, dphi, phi_inv, theta) for a given family and Kendall's tau.
# tau -> theta uses VineCopula::BiCopTau2Par, the same mapping as
# tau_to_theta() in fair-surv.R (used by cif_copula_single), so both agree
# exactly. (Named cge_tau_to_theta to avoid clashing with that function, which
# has the argument order (tau, copula).)
cge_tau_to_theta <- function(copula, tau) {
  if (copula == "clayton") return(BiCopTau2Par(3, tau))
  if (copula == "frank") return(BiCopTau2Par(5, tau))
  stop("unsupported copula: ", copula)
}

arch_gen <- function(copula, tau) {
  if (tau == 0 || copula == "indep") {
    return(list(phi = function(u) -log(u), dphi = function(u) -1 / u,
                phi_inv = function(s) exp(-s), theta = 0))
  }
  th <- cge_tau_to_theta(copula, tau)
  if (copula == "frank") {
    list(phi     = function(u) -log(expm1(-th * u) / expm1(-th)),
         dphi    = function(u) -th / expm1(th * u),
         phi_inv = function(s) -log1p(exp(-s) * expm1(-th)) / th,
         theta = th)
  } else if (copula == "clayton") {
    list(phi     = function(u) (u^(-th) - 1) / th,
         dphi    = function(u) -u^(-th - 1),
         phi_inv = function(s) (1 + th * s)^(-1 / th),
         theta = th)
  }
}

# column-wise cumulative sum of a matrix (fast for n x m with large n)
col_cumsum <- function(M) {
  if (ncol(M) > 1) for (k in 2:ncol(M)) M[, k] <- M[, k - 1] + M[, k]
  M
}

# ---- Closed-form copula-graphic map + linearization -----------------------
# FT, FC : r x m matrices, CIF_T and CIF_C at grid t_1..t_m (t_0 = 0 implicit,
#          CIF(t_0) = 0). Row = unit (Route I) or a single world cell (Route II).
# dT, dC : optional perturbation matrices (same dims as FT after broadcasting).
#          If nrow(FT) == 1 and nrow(dT) = n > 1, FT/FC are broadcast to n rows.
# gen    : output of arch_gen().
# Returns:
#   S    : phi-midpoint estimate S_hat(t_j), r x m (identical to forward recursion)
#   lwr, upr : identified band for S(t_j) (cumulative bound), r x m
#   dS   : linearization grad(g)^T d, same dims as dT (NULL if d not given)
cge_lin <- function(FT, FC, dT = NULL, dC = NULL, gen, eps = 1e-10) {
  FT <- as.matrix(FT); FC <- as.matrix(FC)
  if (!is.null(dT) && nrow(FT) == 1 && nrow(dT) > 1) {
    FT <- FT[rep(1, nrow(dT)), , drop = FALSE]
    FC <- FC[rep(1, nrow(dT)), , drop = FALSE]
  }
  m <- ncol(FT)
  cl <- function(x) pmin(pmax(x, eps), 1)
  FT0 <- cbind(0, FT); FC0 <- cbind(0, FC)
  Stc <- cl(1 - FT[, , drop = FALSE] - FC)                                # S_TC(t_j), j = 1..m
  Stc_prev <- cl(1 - FT0[, -(m + 1), drop = FALSE] - FC0[, -(m + 1), drop = FALSE])  # S_TC(t_{j-1})
  u <- cl(1 - FT - FC0[, -(m + 1), drop = FALSE])      # S_TC(t_{i-1}) - dCIF_{i,T}
  v <- cl(1 - FT0[, -(m + 1), drop = FALSE] - FC)      # S_TC(t_{i-1}) - dCIF_{i,C}
  p <- gen$phi
  Lam <- 0.5 * p(Stc) + 0.5 * col_cumsum(p(u) - p(v))  # phi(S_hat(t_j))
  D   <- p(Stc) - p(u) - p(v) + p(Stc_prev)            # band width increments (>= 0)
  hw  <- 0.5 * col_cumsum(D)                           # half-width on phi-scale
  S   <- gen$phi_inv(Lam)
  out <- list(S = S, lwr = gen$phi_inv(Lam + hw), upr = gen$phi_inv(Lam - hw))
  if (!is.null(dT)) {
    dp <- gen$dphi
    dT0 <- cbind(0, dT); dC0 <- cbind(0, dC)
    term_u <- dp(u) * (dT + dC0[, -(m + 1), drop = FALSE])
    term_v <- dp(v) * (dT0[, -(m + 1), drop = FALSE] + dC)
    dLam <- -0.5 * dp(Stc) * (dT + dC) - 0.5 * col_cumsum(term_u - term_v)
    out$dS <- dLam / dp(S)
  }
  out
}
