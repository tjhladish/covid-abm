
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff-rt.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_rt.png")
))

load(.args[1])

#' comes key'd
rt.dt <- readRDS(.args[2])[
  eval(datefilter)
][, .(
  scenario, realization, date, value = Rt
)]

intscns <- rt.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

plt.dt <- setkeyv(
  rt.dt[scn.dt, on=.(scenario)],
  union(key(rt.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(rt.dt)
gc()

plt.qs <- quantile(
  plt.dt,
  sampleby = "realization",
  probs = qprobs(c(`90`=.9), mid = TRUE, extent = FALSE)
)[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LS", "MS", "HS", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

p <- ggplot(plt.qs) + aes(
  x = date,
  color = act_vac,
  linetype = factor(c("nonpi", "wquar")[quar+1])
) +
  facet_nested(
    cols = vars(talloc), switch = "y", scales = "free_y"
  ) +
  geom_month_background(
    plt.qs, by = "talloc",
    font.size = 3, value.col = "qmed", max.col = "q90h", min.col = "q90l"
  ) +
  geom_ribbon(aes(ymin=q90l, ymax=q90h, fill=act_vac, color=NULL), alpha=0.25) +
  geom_line(aes(y=qmed)) +
  scale_color_strategy() +
  scale_y_continuous(name = "Model Rt") +
  scale_x_null() +
  scale_linetype_quar() +
  scale_alpha(range = c(0.02, 1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    legend.direction = "horizontal"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
