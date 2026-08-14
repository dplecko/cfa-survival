
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
pkgs <- c("data.table", "ggplot2", "randomForestSRC", "zeallot",
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

# rf.cores: OpenMP threads randomForestSRC uses per rfsrc() fit (tree-level parallelism)
# mc.cores: default workers for parallel::mclapply() (fit-level parallelism)
#
# default = many independent fits in parallel (mclapply over N, rf.cores=1 each).
# for a single large fit, set RF_CORES (and usually MC_CORES=1) to instead
# parallelize tree growing within that one fit.
options(
  rf.cores = as.integer(Sys.getenv("RF_CORES", unset = "1")),
  mc.cores = as.integer(Sys.getenv("MC_CORES", unset = n_cores()))
)
data.table::setDTthreads(1)

if (!all(vapply(pkgs, requireNamespace, logical(1L)))) {
  stop("Packages {pkgs} are required in order to proceed.")
  if (!interactive()) q("no", status = 1, runLast = FALSE)
}
for (pkg in pkgs) library(pkg, character.only = TRUE)

