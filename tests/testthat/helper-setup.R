
set.seed(1)
suppressPackageStartupMessages({
  library(data.table)
  library(randomForestSRC)  # for `follic`
  library(survival)
})

# Load project code (non-package)
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
code_dir <- file.path(root, "R")
Rfiles <- list.files(code_dir, full.names = TRUE, pattern = "\\.[Rr]$")
invisible(lapply(Rfiles, sys.source, envir = topenv()))

# Small fixture builder
mk_data <- function(n = 120L) {
  data(follic, package = "randomForestSRC")
  dt <- as.data.table(follic)[1:n]
  dt[, majority := rbinom(.N, 1, 0.5)]
  dt[, rt := NULL]
  dt[, ch := as.integer(ch == "Y")]
  dt[]
}
