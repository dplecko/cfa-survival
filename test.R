
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# select data source
src <- "aics"
out <- "death" # "dcr", "readm"

# prepare the data and the SFM
dat <- load_data("aics", outcome = out)
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")

set.seed(2026)

# for testing
dat_run <- rbind(
  dat[sample(which(dat$majority == 0), size = 500, replace=FALSE)],
  dat[sample(which(dat$majority == 1), size = 500, replace=FALSE)]
)
nboot <- 2

fsurv <- fair_surv(dat_run, X, Z, W, time_var, event_var, 
                   method = "rfs-cf", nboot = nboot,
                   balance_groups = FALSE,
                   split_forest = TRUE,
                   copula = if (out == "readm") "frank" else NULL,
                   tau_grid = if (out == "readm") c(0.1, 0.5, 0.8) else NULL)
