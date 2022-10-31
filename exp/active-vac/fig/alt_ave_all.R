
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_ave_all.png")
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

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(eff.dt)
gc()

plt.qs <- quantile(
  plt.dt,
  j = .(c.averted), sampleby = "realization",
  probs = qprobs(c(`90`=.9, `50`=.5), mid = TRUE, extent = FALSE)
)[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LIC", "MIC", "HIC", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

p <- allplot(
  plt.qs, yl = "Per 10k, Cumulative Averted\nIncidence of ... Relative to Reference",
  withRef = TRUE
)

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
