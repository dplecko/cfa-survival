
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
invisible(lapply(list.files(file.path(root, "r"), full.names = TRUE), source))

res_dir <- file.path(root, "results", "sensitivity")

files <- c(
  classic   = "death_classic.RData",
  irsad_med = "death_irsad_med.RData",
  clin_med  = "death_clin_med.RData"
)

labs <- c(
  classic   = "Current",
  irsad_med = "SES as mediator",
  clin_med  = "Clinical mediators"
)

# load all three results
res <- rbindlist(lapply(names(files), function(z) {
  e <- new.env()
  load(file.path(res_dir, files[[z]]), envir = e)
  x <- as.data.table(e$dr_fsurv$measures)
  x[, zw_split := z]
  x
}))

res[, zw_split := factor(zw_split, levels = names(files), labels = labs)]
res <- res[effect %in% c("ctfde", "ctfie", "ctfse")]

# labels
eff_labs <- c(
  ctfde = "DE",
  ctfie = "IE",
  ctfse = "SE"
)
res[, effect_lab := factor(effect, levels = names(eff_labs), labels = eff_labs)]

# ---------------- markdown table ----------------
# cells: estimate [95% CI]
tab <- copy(res[time_interest %in% c(1, 3, 5, 7, 14, 28, 56, 90, 180)])
tab[, cell := sprintf("%.2f", 100 * value)]

tab <- dcast(
  tab,
  effect_lab + zw_split ~ time_interest,
  value.var = "cell"
)

setnames(tab, c("effect_lab", "zw_split"), c("Ctf-Effect", "Z/W partition"))

knitr::kable(
  tab,
  format = "pipe",
  align = c("l", "l", rep("c", ncol(tab) - 2))
)

# ---------------- plot ----------------
p <- ggplot(
  res[effect != "ctfde"],
  aes(x = time_interest, y = value,
      color = zw_split, fill = zw_split)
) +
  geom_ribbon(
    aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd),
    alpha = 0.15, color = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  facet_wrap(~ effect_lab, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  labs(
    x = "Time",
    y = "Effect",
    color = "Z/W partition",
    fill = "Z/W partition"
  ) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

p

ggsave(
  file.path(res_dir, "zw-sensitivity.png"),
  p, width = 9, height = 4.2
)
