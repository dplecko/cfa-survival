
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), 
                 source))

# select data source
src <- "aics"

# 4-task array: death across all three Z/W splits, readm on the classic split
job_grid <- list(
  list(out = "death",   zw_bnd = "classic"),
  list(out = "death",   zw_bnd = "irsad_med"),
  list(out = "death",   zw_bnd = "clin_med"),
  list(out = "readm", zw_bnd = "classic")
)
task_id <- as.integer(Sys.getenv("SGE_TASK_ID", unset = "1"))
out <- job_grid[[task_id]]$out       # "dcr", "readm"
zw_bnd <- job_grid[[task_id]]$zw_bnd # classic, irsad_med, or clin_med

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

set.seed(2026)

# for testing; runs under SGE (task array) always use the full data
local <- !nzchar(Sys.getenv("SGE_TASK_ID"))
if (local) {

  dat_run <- dat[sample.int(nrow(dat), size = 5 * 10^3, replace=FALSE)]
  nboot <- 10
} else {
  
  dat_run <- dat
  nboot <- 64
}

cat("--- Processing task", task_id, "of", length(job_grid), 
    "for outcome", out, "and Z/W split", zw_bnd, "---\n")

cat("--- Successfully loaded data with", nrow(dat_run), "observations ---\n")

# time grid for the analysis
tgrid <- c(1:10, 14, 28, 56, 90, 180)

# doubly robust estimation
dr_fsurv <- one_step_debias_surv(
  dat_run, X, Z, W, time_var, event_var, time_interest = tgrid,
  copula = if (out == "readm") "frank" else NULL,
  tau_grid = if (out == "readm") c(0.1, 0.5, 0.8),
  # local runs use a 5000-row subsample, so they get their own cache file
  cache_file = f("cache/{out}_{zw_bnd}{if (local) '_local' else ''}.rds"),
  route = if (out == "readm") c("envelope", "II", "I") else "envelope"
)

save(dr_fsurv, file = f("results/sensitivity/{out}_{zw_bnd}.RData"))

autoplot(dr_fsurv)
ggsave(file.path("results", paste0("dr-", out, "-", zw_bnd, ".png")), 
       width = 14, height = 4)

# one plot per sensitivity route (envelope / Route II / Route I) for readm
if (out == "readm") {
  
  for (rt in c("envelope", "II", "I")) {
    
    autoplot(dr_fsurv, route = rt)
    ggsave(file.path("results", paste0("dr-", out, "-", zw_bnd, "-route", rt, 
                                       ".png")), width = 14, height = 4)
  }
}

# library(data.table)
# 
# times_show <- c(7, 28, 90, 180)
# effects_show <- c("tv", "ctfde", "ctfie", "ctfse")
# 
# tab <- copy(dr_fsurv$measures)[
#   event %in% 1:2 & effect %in% effects_show & time_interest %in% times_show
# ]
# 
# tab[, outcome := fifelse(event == 1, "Readmission", "Death")]
# tab[, effect := factor(effect, levels = effects_show,
#                        labels = c("TV", "DE", "IE", "SE"))]
# tab[, value := sprintf("%.2f", 100 * value)]
# 
# tab <- dcast(tab, outcome + effect ~ time_interest, value.var = "value")
# setnames(tab, as.character(times_show), paste0(times_show, " days"))
# 
# knitr::kable(tab, format = "pipe", align = c("l", "l", rep("r", length(times_show))))

# library(data.table)
# library(ggplot2)
# 
# effects_show <- c("tv", "ctfde", "ctfie", "ctfse")
# eff_labs <- c(tv = "TV", ctfde = "Ctf-DE", ctfie = "Ctf-IE", ctfse = "Ctf-SE")
# 
# pd <- copy(dr_fsurv$measures)[event %in% 1:2 & effect %in% effects_show & time_interest <= 180]
# pd[, outcome := factor(fifelse(event == 1, "Readmission", "Death"),
#                        levels = c("Death", "Readmission"))]
# pd[, effect_lab := factor(effect, levels = effects_show, labels = eff_labs)]
# pd[, `:=`(est = value, lo = (value - 1.96 * sd),
#           hi = (value + 1.96 * sd))]
# 
# p <- ggplot(pd, aes(time_interest, est, color = outcome, fill = outcome)) +
#   geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, color = NA) +
#   geom_line(linewidth = .8) +
#   geom_hline(yintercept = 0, linetype = 2, linewidth = .4) +
#   facet_wrap(~effect_lab, nrow = 1, scales = "free_y") +
#   labs(x = "Days after ICU admission", y = "Effect (%)",
#        color = "Outcome", fill = "Outcome") +
#   theme_bw() +
#   scale_y_continuous(labels = scales::percent) +
#   theme(legend.position = "bottom", strip.background = element_blank(),
#         panel.grid.minor = element_blank())
# 
# p
# ggsave(file.path(root, "results", "cr-anzics.png"), p, width = 10, height = 3.2, dpi = 300)
