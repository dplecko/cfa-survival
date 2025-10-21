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
out <- "death" # "dcr", "readm"
balance <- FALSE
splitf <- FALSE

# prepare the data and the SFM
dat <- load_data("aics", outcome = out)
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")

set.seed(2026)

# for testing
local <- FALSE
if (local) {

  dat_run <- dat[sample.int(nrow(dat), size = 1000, replace=FALSE)]
  nboot <- 3
} else {
  
  dat_run <- dat
  nboot <- 64
}

fsurv <- fair_surv(dat_run, X, Z, W, time_var, event_var, 
                   method = "rfs-cf", nboot = nboot,
                   balance_groups = balance,
                   split_forest = splitf,
                   copula = if (out == "readm") "frank" else NULL,
                   tau_grid = if (out == "readm") c(0.1, 0.5, 0.8) else NULL,
                   nodesize = 100)

if (!local) 
  save(fsurv, file = file.path("data", fname(src, out, balance, splitf)))

# load(file.path("data", fname(src, out, balance, splitf)))

# the local analyses need to be adapted!
save <- FALSE
if (local) {

  # (A) Kaplan-Meier Survival Curve Estimates
  km_curves(dat_run)
  if (save) ggsave("results/surv-curves.png", width = 6, height = 4)

  # (B) Survival-TV estimate from fairsurv object
  tv_comparison(fsurv, dat_run)
  if (save) ggsave("results/surv-tv-ot.png", width = 6, height = 4)

  # (C) Survival-TV decomposition with a static comparison
  srv_with_stat(fsurv, dat_run)
  if (save) ggsave("results/decomp-with-static.png", width = 6, height = 4)

  # (D) CHF-TV decomposition with a Cox comparison
  autoplot(fsurv, scale = "chf-ratio") +
    coord_cartesian(ylim=c(0.75, 1.8)) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7),
          legend.box.background = element_rect())

  if (save) ggsave("results/chf-ratio.png", width = 6, height = 4)
}

# load(paste0("data/", src, "_fsurv.RData"))
# load(paste0("data/", src, "_fsurv_", out, ".RData"))
# compare_srv(dat_run, fsurv, "marginal")
# compare_srv(dat, fsurv, "tv")
# ggsave("~/Desktop/split-rf.png", width = 10, height = 4)

