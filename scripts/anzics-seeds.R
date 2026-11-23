root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

# Seed robustness of the ANZICS death analysis: tasks 1 and 3-11 of
# scripts/anzics-survival.R rerun the same analysis on the same data with 10
# seeds (2026, 1, ..., 9), which changes only the cross-fitting folds and the
# random-forest nuisance fits. Compares the spread of the estimates across
# seeds with the per-run standard error, and checks sign/significance stability.
res_dir <- file.path(root, "results", "anzics")
fls <- list.files(res_dir, "^death_classic_seed[0-9]+\\.RData$", full.names = TRUE)
cat(length(fls), "seeds found\n")

res <- rbindlist(lapply(fls, function(fl) {
  e <- new.env()
  load(fl, envir = e)
  x <- as.data.table(e$dr_fsurv$measures)
  x[, seed := as.integer(sub(".*seed([0-9]+)\\.RData$", "\\1", fl))]
  x
}))
res <- res[effect %in% c("tv", "ctfde", "ctfie", "ctfse")]
res[, sig := abs(value) > 1.96 * sd]

# across-seed SD relative to the per-run SE, over all effects x grid days
rat <- res[, .(ratio = sd(value) / median(sd)), by = .(effect, time_interest)]
cat("\n(SD across seeds) / (per-run SE), all effects x days:\n")
print(rat[, .(cells = .N, median = median(ratio), max = max(ratio))], digits = 3)
print(rat[, .(median = median(ratio), max = max(ratio)), by = effect], digits = 3)

# day-7 summary per effect
cat("\nDay 7:\n")
print(res[time_interest == 7, .(
  seeds = .N, mean = mean(value), sd_across = sd(value), per_run_se = median(sd),
  ratio = sd(value) / median(sd), n_pos = sum(value > 0), n_sig = sum(sig),
  z_range = paste(sprintf("%.2f", range(value / sd)), collapse = " to ")
), by = effect], digits = 3)

fwrite(rat, file.path(res_dir, "seed-robustness.csv"))
