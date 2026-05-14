# pkgs
MASS <- requireNamespace("MASS", quietly=TRUE)

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
                     pW = NULL,
                     seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  stopifnot(MASS)
  
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
  
  # T | X,Z,W ~ LogNormal(mu_T, sT)
  muT <- mk_mu(X, Z, W, T_par)
  T <- rlnorm(n, meanlog=muT, sdlog=T_par$sT)
  
  # C | X,Z,W ~ LogNormal(mu_C, sC)
  muC <- mk_mu(X, Z, W, C_par)
  C <- rlnorm(n, meanlog=muC, sdlog=C_par$sT)
  
  M <- pmin(T, C)
  delta <- as.integer(T <= C)
  
  df <- data.frame(majority = X, Z, W, T=T, C=C, event_time=M, event=delta)
  
  list(data=df,
       par=list(T=T_par, C=C_par, A=A, beta=beta, alpha=alpha, B=B,
                SigU=SigU, sZ=sZ, sW=sW, pW=pW),
       cfW = list(Wx0=Wx0, Wx1=Wx1))
}

# ------- survival curves for T -------
S_T_curves <- function(df, par, t_grid) {
  X <- df$X
  Z <- as.matrix(df[ , grepl("^z\\d+$", names(df)), drop=FALSE])
  W <- as.matrix(df[ , grepl("^w\\d+$", names(df)), drop=FALSE])
  muT <- mk_mu(X, Z, W, par$T)
  sT <- par$T$sT
  # n x m matrix of S_T(t)
  M <- outer(muT, t_grid, function(m,t) lnorm_S(t, m, sT))
  dimnames(M) <- list(NULL, paste0("t_", t_grid))
  M
}

# ------- censoring curves for C -------
G_C_curves <- function(df, par, t_grid){
  # allow either df$X or df$majority as treatment column
  X <- if ("X" %in% names(df)) df$X else df$majority
  Z <- as.matrix(df[ , grepl("^z\\d+$", names(df)), drop=FALSE])
  W <- as.matrix(df[ , grepl("^w\\d+$", names(df)), drop=FALSE])
  muC <- mk_mu(X, Z, W, par$C)
  sC <- par$C$sT
  # n x m matrix of G(t) = P(C > t)
  M <- outer(muC, t_grid, function(m,t) lnorm_S(t, m, sC))
  dimnames(M) <- list(NULL, paste0("t_", t_grid))
  M
}

# ground-truth censoring survival G(t | X, Z, W) from generator output
# returns an n x length(t_grid) matrix

ground_truth_G <- function(g, t_grid){
  df <- g$data
  par <- g$par
  G_C_curves(df, par, t_grid)
}

# build a df with new X and W (add)
make_df_with <- function(df, X_new, W_new){
  df2 <- df
  df2$X <- as.integer(X_new)
  wcols <- grepl("^w\\d+$", names(df2))
  df2[, wcols] <- W_new
  df2
}

# potential-outcome survival curves using cfW from generator (add)
S_T_potential_curves_from_gen <- function(g, t_grid){
  df <- g$data; par <- g$par; cfW <- g$cfW
  list(
    S_x0_wx0 = S_T_curves(make_df_with(df, 0L, cfW$Wx0), par, t_grid),
    S_x0_wx1 = S_T_curves(make_df_with(df, 0L, cfW$Wx1), par, t_grid),
    S_x1_wx0 = S_T_curves(make_df_with(df, 1L, cfW$Wx0), par, t_grid),
    S_x1_wx1 = S_T_curves(make_df_with(df, 1L, cfW$Wx1), par, t_grid)
  )
}

ground_truth <- function(g, tg) {
  
  Scf <- S_T_potential_curves_from_gen(g, tg)
  idx1 <- g$data$majority == 1

  tv <- colMeans(Scf$S_x1_wx1[idx1, ]) - colMeans(Scf$S_x0_wx0[!idx1, ])
  
  nde <- colMeans(Scf$S_x1_wx0 - Scf$S_x0_wx0)
  nie <- colMeans(Scf$S_x1_wx0 - Scf$S_x1_wx1)
  nse <- colMeans(Scf$S_x1_wx1[idx1, ]) - colMeans(Scf$S_x1_wx1) -
        (colMeans(Scf$S_x0_wx0[!idx1, ]) - colMeans(Scf$S_x0_wx0))
  
  ctfde <- colMeans( (Scf$S_x1_wx0 - Scf$S_x0_wx0)[!idx1, ] )
  ctfie <- colMeans( (Scf$S_x1_wx0 - Scf$S_x1_wx1)[!idx1, ] )
  ctfse <- colMeans(Scf$S_x1_wx1[!idx1, ]) - colMeans(Scf$S_x1_wx1[idx1, ])
  
  res <- rbind(
    data.frame(time_interest=tg, scale="surv", event=1, effect="tv", value=tv),
    data.frame(time_interest=tg, scale="surv", event=1, effect="nde", value=nde),
    data.frame(time_interest=tg, scale="surv", event=1, effect="nie", value=nie),
    data.frame(time_interest=tg, scale="surv", event=1, effect="nse", value=nse),
    data.frame(time_interest=tg, scale="surv", event=1, effect="ctfde", value=ctfde),
    data.frame(time_interest=tg, scale="surv", event=1, effect="ctfie", value=ctfie),
    data.frame(time_interest=tg, scale="surv", event=1, effect="ctfse", value=ctfse)
  )
  rownames(res) <- NULL
  as.data.table(res)
}

# ------- nested mean: E[ S(x_y, Z, W) | x_w, Z ] via Monte Carlo over W -------

# Generic helper: returns an n x length(t_grid) matrix.
# Monte Carlo draws epsilon_W using the same parametric form as gen_surv.
E_S_xy_given_xw_Z <- function(g, t_grid, xy, xw, B = 200L, seed = NULL){
  if (!is.null(seed)) set.seed(seed)

  df <- g$data
  par <- g$par

  Z <- as.matrix(df[ , grepl("^z\\d+$", names(df)), drop = FALSE])
  q <- length(par$alpha)
  stopifnot(length(par$sW) == q)
  pW <- if (!is.null(par$pW)) par$pW else runif(1, 0.2, 0.4)
  stopifnot(length(pW) == 1)

  n <- nrow(df)
  m <- length(t_grid)
  out <- matrix(0, n, m)

  muW_xw <- outer(rep.int(xw, n), par$alpha) + Z %*% par$B

  for (b in seq_len(B)) {
    Eb <- matrix(rnorm(n*q, 0, rep(par$sW, each = n)), n, q)
    Mb <- rbinom(n, 1, pW)
    Wb <- (Z %*% par$B) + Eb + outer(as.integer(xw) * Mb, par$alpha)
    Sb <- S_T_curves(make_df_with(df, as.integer(xy), Wb), par, t_grid)
    out <- out + Sb
  }

  out / B
}

# Convenience wrappers (four combinations)
E_S_x0_given_x0_Z <- function(g, t_grid, B = 200L, seed = NULL)
  E_S_xy_given_xw_Z(g, t_grid, xy = 0L, xw = 0L, B = B, seed = seed)

E_S_x0_given_x1_Z <- function(g, t_grid, B = 200L, seed = NULL)
  E_S_xy_given_xw_Z(g, t_grid, xy = 0L, xw = 1L, B = B, seed = seed)

E_S_x1_given_x0_Z <- function(g, t_grid, B = 200L, seed = NULL)
  E_S_xy_given_xw_Z(g, t_grid, xy = 1L, xw = 0L, B = B, seed = seed)

E_S_x1_given_x1_Z <- function(g, t_grid, B = 200L, seed = NULL)
  E_S_xy_given_xw_Z(g, t_grid, xy = 1L, xw = 1L, B = B, seed = seed)
