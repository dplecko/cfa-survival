#!/u/local/apps/R/4.2.2/gcc-4.8.5_intel-2020.4/bin/Rscript
#$ -cwd
#$ -j y
#$ -N fairsurv
#$ -pe shared 8
#$ -l h_rt=1:00:00,h_data=1G

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
local <- FALSE
if (local) {

  dat_run <- dat[sample.int(nrow(dat), size = 10^3, replace=FALSE)]
  nboot <- 10
} else {
  
  dat_run <- dat
  nboot <- 64
}

# time grid for the analysis
tgrid <- c(1:10, 14, 28, 56, 90, 180)

# model-based estimation
fsurv <- fair_surv(dat_run, X, Z, W, time_var, event_var,
                   method = "rfs-cf", nboot = 1,
                   copula = if (out == "readm") "frank" else NULL,
                   time_interest = 100,
                   tau_grid = if (out == "readm") c(0.1, 0.5, 0.8) else NULL)


# doubly robust estimation
dr_fsurv <- one_step_debias_surv(
  dat_run, X, Z, W, time_var, event_var, time_interest = tgrid,
  copula = if (out == "readm") "frank" else NULL,
  tau_grid = if (out == "readm") c(0.1, 0.5, 0.8)
)


autoplot(dr_fsurv)
ggsave(file.path("results", paste0("dr-", out, ".png")), width = 14, height = 4)

