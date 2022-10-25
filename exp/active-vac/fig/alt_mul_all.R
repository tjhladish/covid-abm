
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_eff_all.png")
))

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  eval(datefilter) & eval(outfilter)
][, .(
  scenario, realization, date, outcome, c.effectiveness
)]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(eff.dt)
gc()

plt.dt[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LIC", "MIC", "HIC", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

# eff = 1 - mult => mult = 1-eff

# inverting here, because want <1 on upper plane of plot
plt.dt[, c.mult := 1-c.effectiveness ]
plt.dt[, l.c.mult := -log(c.mult, 2) ]

plt.qs <- quantile(
  plt.dt,
  j = .(l.c.mult), sampleby = "realization",
  probs = qprobs(c(`90`=.9, `50`=.5), mid = TRUE, extent = FALSE)
)[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LIC", "MIC", "HIC", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

muls <- c(1, 5/4, 4/3, 3/2)
muls <- log(unique(c(rev(1/muls), muls)), 2)

p <- ggplot(plt.qs) + aes(
  x = date,
  color = act_vac,
  linetype = factor(c("nonpi", "wquar")[quar+1])
) +
  facet_nested(
    rows = vars(outcome), cols = vars(talloc), switch = "y",
    labeller = labeller(
      outcome = c(inf = "Infection", sev = "Severe Disease", deaths = "Deaths")
    )
  ) +
  geom_month_background(
    plt.qs, by = c(row="outcome", col="talloc"),
    font.size = 3, value.col = "qmed", max.col = "q90h", min.col = "q90l",
    ymax = .6, ymin = -.6, m.y = 0.95, y.y = 0.91
  ) +
  geom_ribbon(aes(ymin=q90l, ymax=q90h, fill=act_vac, color=NULL), alpha=0.15) +
  #  geom_ribbon(aes(ymin=q50l, ymax=q50h, fill=act_vac, color=NULL), alpha=0.25) +
  geom_line(aes(y=qmed)) +
  scale_color_strategy() +
  scale_y_continuous(
    name = "Cumulative Relative\nMultiplier in Incidence of ... (log scale)",
    breaks = muls,
    labels = c("3/2x", "4/3x", "5/4x", "1x", "4/5x", "3/4x", "2/3x")
  ) +
  coord_cartesian(ylim = c(-.6, .6), expand = FALSE) +
  scale_x_null() +
  scale_linetype_quar() +
  #  scale_alpha(range = c(0.02, 1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    legend.direction = "horizontal"
  ) +
  geom_hline(
    aes(yintercept = 0, color = "none", linetype = "nonpi"),
    show.legend = FALSE, data = \(dt) dt[,.SD[1],by=.(outcome, talloc)]
  ) +
  geom_texthline(
    aes(yintercept = 0, color = "none", label = "Reference\nProgram"),
    inherit.aes = FALSE, show.legend = FALSE, data = \(dt) dt[talloc == "LIC",.SD[1],by=.(outcome, talloc)],
    hjust = 0, gap = FALSE
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
