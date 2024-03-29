
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds", "vocwindows.rds")),
  file.path("fig", "output", "alt_ave_all_ns.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  eval(datefilter) & eval(outfilter)
][, .(
  scenario, realization, date, outcome, averted
)]

eff.dt[order(date), c.averted := cumsum(averted), by=.(scenario, realization, outcome)]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns][eval(seasfilter)]

takeover.wins <- readRDS(.args[4])

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario), nomatch = 0],
  union(key(eff.dt), colnames(scn.dt))
)[inf_con == FALSE][eval(seasfilter)]

rm(eff.dt)
gc()

plt.qs <- plt.prep(plt.dt, j = .(c.averted))

mm.ref <- plt.qs[,.(ymin = min(q90l), ymax = max(q90h)),by=.(outcome)]
mm.ref[, ymin := ymin - .15*(ymax-ymin) ]
mm.ref[, yspan := ymax - ymin ]
tw <- takeover.wins[q == 0.5][CJ(measure, outcome = mm.ref$outcome), on=.(measure)][mm.ref, on=.(outcome)]
tw[, end := pmin(end, plt.qs[, max(date)])]
tw[, mid := start + (end-start)/2 ]

p <- allplot(
  plt.qs, yl = "Per 10k, Cumulative Averted\nIncidence of ... Relative to Reference",
  withRef = TRUE, ins = voc.band(tw)
)

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
