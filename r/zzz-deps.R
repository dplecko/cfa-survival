
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
pkgs <- c("data.table", "ggplot2", "ricu", "randomForestSRC", "zeallot",
          "xgboost", "survival", "parallel", "abind", "VineCopula",
          "assertthat")

Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))
Sys.setenv("RICU_SRC_LOAD" = 
             "mimic,miiv,aumc,hirid,eicu,eicu_demo,mimic_demo,anzics,sic")

n_cores <- function() {
  
  as.integer(
    Sys.getenv("NSLOTS", unset = parallel::detectCores() - 1)
  )
}

options(rf.cores = 1, mc.cores = n_cores())   # rfsrc=1 inner; mclapply=N outer
data.table::setDTthreads(1)

if (!all(vapply(pkgs, requireNamespace, logical(1L)))) {
  stop("Packages {pkgs} are required in order to proceed.")
  if (!interactive()) q("no", status = 1, runLast = FALSE)
}
for (pkg in pkgs) library(pkg, character.only = TRUE)

