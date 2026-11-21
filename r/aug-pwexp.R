# Censoring augmentation under piecewise-exponential interpolation.
#
# The nuisances S (all-cause survival), G (censoring survival) and CIF_j are
# only available on tgrid. Within each cell (t_{k-1}, t_k] all hazards are
# taken to be constant, and this single interpolation is used everywhere a
# nuisance is needed off the grid (jump terms at M, IPCW weight G(M)), while
# the compensator integrals are computed in closed form under it. The
# jump term and its compensator then cancel in expectation exactly for the
# interpolated nuisances (Lemmas 1/2), with no quadrature.

AUG_VERSION <- "pwexp-v1"

# lower clip applied once to S and G, before anything else is derived
PW_EPS <- 1e-6

# cell geometry
pw_cells <- function(tgrid) list(t0 = c(0, head(tgrid, -1)), h = diff(c(0, tgrid)))

# n x m matrix of cell hazard rates from survival-type values V (n x m, on tgrid)
pw_rate <- function(V, h) {
  Vp <- cbind(1, V[, -ncol(V), drop = FALSE])
  r  <- pmin(pmax(V / pmax(Vp, 1e-12), 1e-12), 1)
  sweep(-log(r), 2, h, "/")
}

# cell index of M: k such that t_{k-1} < M <= t_k  (clamped to 1..m)
pw_cell <- function(M, tgrid) {
  k <- findInterval(M, c(0, tgrid), left.open = TRUE)
  pmin(pmax(k, 1L), length(tgrid))
}

# value of V (n x m, V(0) = 1) at M under piecewise-exponential interpolation
# (for M > t_m the cell offset is capped at h_m, i.e. V(t_m) is returned)
pw_at <- function(V, rate, M, tgrid, cells) {
  k  <- pw_cell(M, tgrid); i <- seq_along(M)
  Vp <- ifelse(k == 1L, 1, V[cbind(i, pmax(k - 1L, 1L))])
  d  <- pmin(pmax(M - cells$t0[k], 0), cells$h[k])
  Vp * exp(-rate[cbind(i, k)] * d)
}

# per-cell at-risk integrals (n x m):
#   IGS_ik = int_{t0_k}^{t0_k + D_ik} du / (S(u) G(u)),  IG_ik = int ... du / G(u)
# with D_ik = min(max(M_i - t0_k, 0), h_k)
pw_integrals <- function(S, G, a, b, M, cells) {
  n <- length(M); m <- length(cells$h)
  Sp <- cbind(1, S[, -m, drop = FALSE]); Gp <- cbind(1, G[, -m, drop = FALSE])
  D  <- pmin(pmax(outer(M, cells$t0, "-"), 0), matrix(cells$h, n, m, byrow = TRUE))
  cc <- a + b
  fGS <- ifelse(cc < 1e-10, D, expm1(cc * D) / pmax(cc, 1e-300))
  fG  <- ifelse(a  < 1e-10, D, expm1(a  * D) / pmax(a,  1e-300))
  list(IGS = fGS / (Sp * Gp), IG = fG / Gp)
}

# row-wise cumulative sum of an n x m matrix (column loop; fast for large n)
pw_rowcumsum <- function(X) {
  if (ncol(X) > 1L) for (k in 2:ncol(X)) X[, k] <- X[, k] + X[, k - 1L]
  X
}

# Pseudo-outcomes on tgrid for one set of test observations.
#   S:     n x m survival (NIC) or all-cause survival (CR) on tgrid
#   G:     n x m censoring survival on tgrid
#   cif:   n x m x J array of CIFs (CR), NULL for NIC
#   chf:   n x m censoring CHF; unused here (kept so that the pre-fix
#          implementation can share this signature)
#   M, delta: observed times and event indicators (0 = censored)
# Returns list(H, clip): H is n x m (NIC) or n x m x J (CR); clip counts the
# entries affected by the numerical guards.
pw_augment <- function(S, G, cif, chf, M, delta, tgrid, martingale_debias) {

  is_cr <- !is.null(cif)
  n <- length(M); m <- length(tgrid)
  S <- matrix(S, n, m); G <- matrix(G, n, m)

  clip <- list(
    S_floor = sum(S < PW_EPS), G_floor = sum(G < PW_EPS),
    S_nonmono = 0L, G_nonmono = 0L, pi_zero = 0L, n_entries = n * m
  )

  # clip once, at the source: everything below is derived from these
  S <- pmax(S, PW_EPS); G <- pmax(G, PW_EPS)
  clip$S_nonmono <- sum(S > cbind(1, S[, -m, drop = FALSE]))
  clip$G_nonmono <- sum(G > cbind(1, G[, -m, drop = FALSE]))

  cells <- pw_cells(tgrid)
  a   <- pw_rate(G, cells$h)
  G_M <- pw_at(G, a, M, tgrid, cells)

  atrisk_after <- outer(M, tgrid, ">")    # 1(M > t_k)
  obs_by       <- !atrisk_after           # 1(M <= t_k)

  if (martingale_debias) {
    b    <- pw_rate(S, cells$h)
    S_M  <- pw_at(S, b, M, tgrid, cells)
    I    <- pw_integrals(S, G, a, b, M, cells)
    cumA <- pw_rowcumsum(a * I$IGS)
    cens_by <- obs_by & (delta == 0)      # 1(M <= t_k, delta = 0)
  }

  if (!is_cr) {

    H <- atrisk_after / G
    if (martingale_debias)
      H <- H + S * cens_by / (S_M * G_M) - S * cumA
  } else {

    J <- dim(cif)[3]
    cif <- array(cif, c(n, m, J))
    H <- array(NA_real_, c(n, m, J))
    Sp <- cbind(1, S[, -m, drop = FALSE])
    if (martingale_debias) {
      kM  <- pw_cell(M, tgrid)
      idx <- cbind(seq_len(n), kM)
      dS  <- Sp - S
    }

    for (j in seq_len(J)) {

      H[, , j] <- (obs_by & (delta == j)) / G_M

      if (martingale_debias) {
        Fj  <- cif[, , j]
        Fjp <- cbind(0, Fj[, -m, drop = FALSE])
        ok  <- dS > 1e-12
        pij <- ifelse(ok, (Fj - Fjp) / ifelse(ok, dS, 1), 0)
        clip$pi_zero <- clip$pi_zero + sum(!ok)

        Bj   <- pw_rowcumsum(a * ((Fjp + pij * Sp) * I$IGS - pij * I$IG))
        Fj_M <- Fjp[idx] + pij[idx] * (Sp[idx] - S_M)

        zeta1 <- cens_by / G_M * (Fj - Fj_M) / S_M
        zeta2 <- Fj * cumA - Bj
        H[, , j] <- H[, , j] + zeta1 - zeta2
      }
    }
  }

  list(H = H, clip = clip)
}
