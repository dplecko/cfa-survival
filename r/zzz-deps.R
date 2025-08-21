
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
pkgs <- c("data.table", "ggplot2", "ricu", "randomForestSRC", "zeallot",
          "xgboost", "ranger", "survival", "parallel", "abind", "VineCopula")

Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))
Sys.setenv("RICU_SRC_LOAD" = 
             "mimic,miiv,aumc,hirid,eicu,eicu_demo,mimic_demo,anzics,sic")

n_cores <- function() {
  
  as.integer(
    Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()-1)
  )
}

options(rf.cores=n_cores(), mc.cores=n_cores(), ranger.num.threads=n_cores())

if (!all(vapply(pkgs, requireNamespace, logical(1L)))) {
  stop("Packages {pkgs} are required in order to proceed.")
  if (!interactive()) q("no", status = 1, runLast = FALSE)
}
for (pkg in pkgs) library(pkg, character.only = TRUE)

