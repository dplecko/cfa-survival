
file_report <- function(fl) {
  
  agg <- readRDS(file.path(root, "results", fl))$agg
  agg_po <- readRDS(file.path(root, "results", fl))$agg_po
  
  cat(nrow(agg) / (100 * 2 * 4), "seeds in the run\n")
  
  psi001_true <- agg_po[time_interest == eval_grid[2] & xz == 0 & xw == 0 & xy == 1 & event == 2]$truth
  psi011_hat <- agg_po[time_interest == eval_grid[2] & xz == 0 & xw == 1 & xy == 1 & event == 2]$value
  psi001_hat <- agg_po[time_interest == eval_grid[2] & xz == 0 & xw == 0 & xy == 1 & event == 2]$value
  psi011_true <- agg_po[time_interest == eval_grid[2] & xz == 0 & xw == 1 & xy == 1 & event == 2]$truth
  
  
  # format this part
  sd011   <- sd(psi011_hat)
  sd001   <- sd(psi001_hat)
  bias011 <- mean(psi011_hat - psi011_true)
  bias001 <- mean(psi001_hat - psi001_true)
  
  ie_hat  <- psi001_hat - psi011_hat
  ie_true <- psi001_true - psi011_true
  sd_ie   <- sd(ie_hat)
  bias_ie <- mean(ie_hat - ie_true)
  
  cat(sprintf(
    paste0(
      "psi011: SD = %.6f | bias = %+.6f | bias/SD = %+.3f sigma\n",
      "psi001: SD = %.6f | bias = %+.6f | bias/SD = %+.3f sigma\n",
      "corr(psi001, psi011) = %.4f\n",
      "Ctf-IE: SD = %.6f | bias = %+.6f | bias/SD = %+.3f sigma\n"
    ),
    sd011, bias011, bias011 / sd011,
    sd001, bias001, bias001 / sd001,
    cor(psi001_hat, psi011_hat),
    sd_ie, bias_ie, bias_ie / sd_ie
  ))
  
  
  s_po <- readRDS(file.path(root, "results", fl))$summary_po
  tab <- dcast(s_po[time_interest %in% eval_grid],
               event + xz + xw + xy ~ time_interest, value.var = "coverage")
  for (j in names(tab)[-(1:4)]) {
    tab[[j]] <- sprintf("%.0f%%", 100 * tab[[j]])
  }
  print(knitr::kable(tab, format = "pipe",
                     align = c("r", "l", rep("r", length(eval_grid)))))
  
  s <- readRDS(file.path(root, "results", fl))$summary
  tab <- dcast(s[time_interest %in% eval_grid],
               event + effect ~ time_interest, value.var = "coverage")
  for (j in names(tab)[-(1:2)]) {
    tab[[j]] <- sprintf("%.0f%%", 100 * tab[[j]])
  }
  print(knitr::kable(tab, format = "pipe",
                     align = c("r", "l", rep("r", length(eval_grid)))))
  
  invisible(TRUE)
}

fl <- "cr-dml-weibull-10000-zeta-96s.rds"
file_report(fl)


fl <- file.path(root, "results", fl)
decomp <- readRDS(fl)$decomp
agg_po <- readRDS(fl)$agg_po
agg <- readRDS(fl)$agg

pot <- 
  agg_po[xz == 0 & xw == 0 & xy == 1 & event == 2 & time_interest == eval_grid[2]]




qqplot(pot$A0 / sd(pot$A0), rnorm(10^4))
abline(0, 1)

# this confirms, A0 is N(0, var(phi))

quantile(pot$A0 / sd(pot$A0), seq(0, 1, 0.05))

quantile(pot$A1 / sd(pot$A0))


quantile(pot$e2_t2 / sd(pot$A0))

quantile(pot$e2_t3 / sd(pot$A0))

agg_po[xz == 0 & xw == 0 & xy == 1 & event == 2 & time_interest == eval_grid[2]]$e2_t2 -
  agg_po[xz == 0 & xw == 1 & xy == 1 & event == 2 & time_interest == eval_grid[2]]$e2_t2


mean(
  agg_po[xz == 0 & xw == 0 & xy == 1 & event == 2 & time_interest == eval_grid[2]]$e2_t3 -
    agg_po[xz == 0 & xw == 1 & xy == 1 & event == 2 & time_interest == eval_grid[2]]$e2_t3
) / 
  sd(agg[effect == "ctfie" & event == 2 & time_interest == eval_grid[2]]$value)

