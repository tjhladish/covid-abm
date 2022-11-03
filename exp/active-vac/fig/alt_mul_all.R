
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_mul_all.png")
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

p <- allplot(
  plt.qs, yl = "Cumulative Relative\nOutcome Multiplier of ...",
  withRef = TRUE
) + coord_cartesian(ylim = c(-0.75, 0.75), clip = "off") +
  scale_y_continuous(
    name = "Cumulative Relative\nMultiple of ...", breaks = (-3:3)/4,
    labels = \(x) sprintf("%.1fx",(2^-x))
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
