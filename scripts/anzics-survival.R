
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# select data source
src <- "aics"

# 13-task array (run as a Slurm array, e.g. sbatch --array=1-13):
#   1-2: main results (death, readm); 3-11: death seed robustness (task 1 is
#   the 10th seed; summarized by scripts/anzics-seeds.R); 12-13: Z/W split
#   sensitivity for death (figure: scripts/anzics-zw-sensitivity.R)
job_grid <- c(
  list(list(out = "death", zw_bnd = "classic", seed = 2026),
       list(out = "readm", zw_bnd = "classic", seed = 2026)),
  lapply(1:9, function(s) list(out = "death", zw_bnd = "classic", seed = s)),
  list(list(out = "death", zw_bnd = "irsad_med", seed = 2026),
       list(out = "death", zw_bnd = "clin_med", seed = 2026))
)
# Slurm array id (falls back to Grid Engine); no id = local run on a subsample
task_env <- Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("SGE_TASK_ID", ""))
task_id  <- if (nzchar(task_env)) as.integer(task_env) else 1L
local    <- !nzchar(task_env) || nzchar(Sys.getenv("ANZICS_LOCAL"))  # ANZICS_LOCAL=1: 5k subsample
job <- job_grid[[task_id]]
out <- job$out       # "death", "readm"
zw_bnd <- job$zw_bnd # classic, irsad_med, or clin_med
seed <- job$seed

# prepare the data and the SFM
dat <- load_data("aics", outcome = out)
c(X, Z, W, event_var, time_var) %<-% attr(dat, "sfm")
if (zw_bnd == "irsad_med") {
  
  Z <- c("age", "sex")
  W <- c("irsad", "anz_cmb", "frailty", "apache_iii_diag", "elective", 
         "apache_iii_rod")
} else if (zw_bnd == "clin_med") {
  
  Z <- c("age", "sex", "irsad", "anz_cmb", "frailty")
  W <- c("apache_iii_diag", "elective", "apache_iii_rod")
}

set.seed(seed)

# for testing; runs under Slurm/SGE (task array) always use the full data
if (local) {

  dat_run <- dat[sample.int(nrow(dat), size = 5 * 10^3, replace=FALSE)]
  nboot <- 10
} else {
  
  dat_run <- dat
  nboot <- 64
}

cat("--- Processing task", task_id, "of", length(job_grid), 
    "for outcome", out, "and Z/W split", zw_bnd, "seed", seed, "---\n")

cat("--- Successfully loaded data with", nrow(dat_run), "observations ---\n")

# time grid for the analysis
tgrid <- c(1:10, 14, 21, 28, 42, 56, 70, 90, 120, 150, 180)

# doubly robust estimation
dr_fsurv <- one_step_debias_surv(
  dat_run, X, Z, W, time_var, event_var, time_interest = tgrid,
  copula = if (out == "readm") "frank" else NULL,
  tau_grid = if (out == "readm") c(0.1, 0.5, 0.8),
  # one cache per job (seed in the name: seeds must not share a cross-fit);
  # local runs use a 5000-row subsample, so they get their own cache file
  cache_file = f("cache/anzics/{out}_{zw_bnd}_seed{seed}{if (local) '_local' else ''}.rds"),
  route = if (out == "readm") c("envelope", "II", "I") else "envelope"
)

res_dir <- file.path(root, "results", "anzics")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
tag <- f("{out}_{zw_bnd}_seed{seed}{if (local) '_local' else ''}")
save(dr_fsurv, file = file.path(res_dir, f("{tag}.RData")))

# plots only for the primary seed (seed-robustness runs are tables only);
# wrapped in try() so a plotting failure never fails the job (results saved)
if (seed == 2026) try({

  autoplot(dr_fsurv)
  ggsave(file.path(res_dir, f("dr-{tag}.png")), width = 14, height = 4)

  # one plot per sensitivity route (envelope / Route II / Route I) for readm
  if (out == "readm") {

    for (rt in c("envelope", "II", "I")) {

      autoplot(dr_fsurv, route = rt)
      ggsave(file.path(res_dir, f("dr-{tag}-route{rt}.png")), width = 14, height = 4)
    }
  }
})

# competing-risks figure (readmission and death as competing events; Fig. in
# the ANZICS appendix): classical measures for both causes over 180 days
if (out == "readm" && seed == 2026) try({

  effects_show <- c("tv", "ctfde", "ctfie", "ctfse")
  eff_labs <- c(tv = "TV", ctfde = "Ctf-DE", ctfie = "Ctf-IE", ctfse = "Ctf-SE")

  pd <- copy(as.data.table(dr_fsurv$measures))[
    event %in% 1:2 & effect %in% effects_show & time_interest <= 180
  ]
  pd[, outcome := factor(fifelse(event == 1, "Readmission", "Death"),
                         levels = c("Death", "Readmission"))]
  pd[, effect_lab := factor(effect, levels = effects_show, labels = eff_labs)]
  pd[, `:=`(lo = value - 1.96 * sd, hi = value + 1.96 * sd)]

  p <- ggplot(pd, aes(time_interest, value, color = outcome, fill = outcome)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, color = NA) +
    geom_line(linewidth = .8) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = .4) +
    facet_wrap(~effect_lab, nrow = 1, scales = "free_y") +
    labs(x = "Days after ICU admission", y = "Effect (%)",
         color = "Outcome", fill = "Outcome") +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "bottom", strip.background = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(file.path(res_dir, f("cr-anzics{if (local) '_local' else ''}.png")),
         p, width = 10, height = 3.2, dpi = 300)
})
