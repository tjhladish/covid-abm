
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest.rds", "digest-key.rds", "vocwindows.rds")),
  file.path("fig", "output", "alt_inc_all_ns.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

intfilter <- expression(realization >= 0)

#' comes key'd
inc.dt <- readRDS(.args[2])[
  eval(datefilter) & eval(intfilter) & eval(outfilter)
][, .(
  scenario, realization, date, outcome, value
)]

intscns <- inc.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[eval(seasfilter)]

takeover.wins <- readRDS(.args[4])

intscn.dt <- scn.dt[scenario %in% intscns]
# reconstructing reference scenarios
refscn.dt <- scn.dt[quar == FALSE & pas_vac == TRUE & act_vac == "none"]

# incref.dt <- inc.dt[
#   intscn.dt, on=.(scenario)
# ][(act_vac == "ring") & (quar == FALSE)][ # only need to go from one reference
#   refscn.dt, on =.(act_alloc = pas_alloc, inf_con)
# ][,.(
#   value = value[1]
# ), by=.(scenario = i.scenario, realization, date, outcome)
# ][refscn.dt, on=.(scenario)]

plt.dt <- setkeyv(
  inc.dt[intscn.dt, on=.(scenario)],
  union(key(inc.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(inc.dt)
gc()

plt.qs <- plt.prep(plt.dt, j = .(value))

mm.ref <- plt.qs[,.(ymin = min(q90l), ymax = max(q90h)),by=.(outcome)]
mm.ref[, ymin := ymin - .15*(ymax-ymin) ]
mm.ref[, yspan := ymax - ymin ]
tw <- takeover.wins[q == 0.5][CJ(measure, outcome = mm.ref$outcome), on=.(measure)][mm.ref, on=.(outcome)]
tw[, end := pmin(end, plt.qs[, max(date)])]
tw[, mid := start + (end-start)/2 ]

p <- allplot(
  plt.qs, yl = "Per 10k, Incidence of ...",
  withRef = FALSE, ins = voc.band(tw)
)

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
