#!/burg/opt/R-4.3.1/bin/Rscript
#SBATCH --job-name=fairsurv
#SBATCH --account=dsi
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --time=48:00:00

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# select data source
src <- "aics"
out <- "readm" # "dcr", "readm"

# prepare the data and the SFM
dat <- load_data("aics", outcome = out)
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")

set.seed(2026)

# for testing
local <- TRUE
if (local) {

  dat_run <- dat[sample.int(nrow(dat), size = 10^3, replace=FALSE)]
  nboot <- 10
} else {
  
  dat_run <- dat
  nboot <- 64
}

tgrid <- c(1:10, 14, 28, 56, 90, 180)

fsurv <- fair_surv(dat_run, X, Z, W, time_var, event_var, 
                   method = "rfs-cf", nboot = 1,
                   copula = if (out == "readm") "frank" else NULL,
                   time_interest = 100,
                   tau_grid = if (out == "readm") c(0.1, 0.5, 0.8) else NULL)


t0 <- Sys.time()
dr_est <- one_step_debias_surv(
  dat_run, X, Z, W, time_var, event_var, time_interest = tgrid,
  copula = if (out == "readm") "frank" else NULL,
  tau_grid = if (out == "readm") c(0.1, 0.5, 0.8)
)
Sys.time() - t0


ggplot(dr_est, aes(x = time_interest, y = value, color = effect, fill = effect)) +
  geom_line() + theme_bw() +
  geom_ribbon(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd),
              alpha = 0.4, linewidth = 0) +
  facet_wrap(~effect, scales = "free")


if (!local) 
  save(fsurv, file = file.path("data", fname(src, out, balance, splitf)))

# load(file.path("data", fname(src, out, balance, splitf)))

# the local analyses need to be adapted!
save <- FALSE
if (local) {

  # (C) Survival-TV decomposition with a static comparison
  srv_with_stat(fsurv, dat_run)
  if (save) ggsave("results/decomp-with-static.png", width = 6, height = 4)

  
}

