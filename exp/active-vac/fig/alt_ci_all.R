
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds", "vocwindows.rds")),
  file.path("fig", "output", "alt_ci_all_ns.png")
))

load(.args[1])

intfilter <- if (interactive()) expression(realization < 10) else expression(realization >= 0)

#' comes key'd
inc.dt <- readRDS(.args[2])[
  eval(intfilter) & eval(outfilter)
][
  eval(datefilter)
][, .(
  scenario, realization, date, outcome, value, averted
)]

inc.dt[order(date),
  c("c.value", "c.averted") := .(
    cumsum(value), cumsum(averted)
  ), by=.(scenario, realization, outcome)
][, i.c.value := c.value + c.averted ]

intscns <- inc.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[eval(seasfilter)]

takeover.wins <- readRDS(.args[4])

intscn.dt <- scn.dt[scenario %in% intscns]
# reconstructing reference scenarios
refscn.dt <- scn.dt[quar == FALSE & pas_vac == TRUE & act_vac == "none"]

incref.dt <- copy(inc.dt)[
  intscn.dt, on=.(scenario)
][(act_vac == "ring") & (quar == FALSE)][ # only need to go from one reference
  refscn.dt, on =.(act_alloc = pas_alloc, inf_con, quar)
][,.(
  c.value = i.c.value
), by=.(i.scenario, realization, date, outcome)
][refscn.dt, on=.(i.scenario = scenario)]
setnames(incref.dt, "i.scenario", "scenario")

plt.dt <- setkeyv(
  rbind(
    inc.dt[, c(key(inc.dt), "c.value"), with = FALSE][intscn.dt, on=.(scenario)],
    incref.dt
  ),
  union(key(inc.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(inc.dt)
gc()

plt.qs <- plt.prep(plt.dt, j = .(c.value))

mm.ref <- plt.qs[,.(ymin = min(q90l), ymax = max(q90h)),by=.(outcome)]
mm.ref[, ymin := ymin - .15*(ymax-ymin) ]
mm.ref[, yspan := ymax - ymin ]
tw <- takeover.wins[q == 0.5][CJ(measure, outcome = mm.ref$outcome), on=.(measure)][mm.ref, on=.(outcome)]
tw[, end := pmin(end, plt.qs[, max(date)])]
tw[, mid := start + (end-start)/2 ]

p <- allplot(
  plt.qs, yl = "Per 10k, Cumulative\nIncidence of ...",
  withRef = FALSE, ins = voc.band(tw)
)

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
