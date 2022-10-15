
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_ci_all.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#' comes key'd
inc.dt <- readRDS(.args[2])[
  eval(datefilter)
][
  outcome %in% c("inf", "sev", "deaths")
][, .(
  scenario, realization, date, outcome, value, averted
)]

inc.dt[order(date),
  c("c.value", "i.c.value") := .(
    cumsum(value), cumsum(value + averted)
  ), by=.(scenario, realization, outcome)
]

intscns <- inc.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])

intscn.dt <- scn.dt[scenario %in% intscns]
# reconstructing reference scenarios
refscn.dt <- scn.dt[quar == FALSE & pas_vac == TRUE & act_vac == "none"]

incref.dt <- inc.dt[
  intscn.dt, on=.(scenario)
][
  refscn.dt, on =.(act_alloc = pas_alloc, inf_con)
][,.(
  c.value = i.c.value[1]
), by=.(scenario = i.scenario, realization, date, outcome)
][refscn.dt, on=.(scenario)]


plt.dt <- setkeyv(
  rbind(
    inc.dt[, .SD, .SDcols = -c("value", "i.c.value", "averted")][intscn.dt, on=.(scenario)],
    incref.dt
  ),
  union(key(inc.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(inc.dt)
gc()

plt.qs <- quantile(
  plt.dt,
  j = .(c.value), sampleby = "realization",
  probs = qprobs(c(`90`=.9), mid = TRUE, extent = FALSE)
)[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LIC", "MIC", "HIC", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

p <- allplot(
  plt.qs, ylab = "Per 10k, Cumulative\nIncidence of ...",
  withRef = FALSE
)

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
