
.pkgs <- c(
  "data.table", "cabputils"
)

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "process", "preference.rds")
))

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  eval(datefilter) & eval(outfilter)
][, .(
  scenario, realization, date, outcome, value, averted
)]

eff.dt[order(date), cvalue := cumsum(value), by = .(scenario, realization, outcome)]
# need to reconstruct reference scenarios
# eff.dt has no information on "do nothing" scenarios
# instead, it is baselined against passive vaccination, without quarantine,
# for the relevant supply levels + w/ vs w/o conditioning vaccination on infection
# these can be reconstructed from any of their matches + averted column

scn.dt <- readRDS(.args[3])
act_levels <- scn.dt[, levels(act_vac)]
act_levels <- c(act_levels[1], "passive", act_levels[-1])
scn.dt[, act_vac := factor(
  fifelse(pas_vac, "passive", as.character(act_vac))
, levels = act_levels, ordered = TRUE)]

nons <- setdiff(1:max(scn.dt$scenario), eff.dt[, unique(scenario)])
refs <- scn.dt[scenario %in% nons & pas_vac == TRUE]
refs[, pas_vac := FALSE ][, act_vac := "ring"][, act_alloc := pas_alloc][, pas_alloc := "none"][, oscen := scenario]
refs$scenario <- NULL
refs[scn.dt, on = .(quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con, season), scenario := scenario]

refeff.dt <- eff.dt[scenario %in% refs$scenario]
refeff.dt[, value := value + averted][, averted := 0]
refeff.dt[order(date), cvalue := cumsum(value), by = .(scenario, realization, outcome)]
refeff.dt[refs[, .(scenario, oscen)], on = .(scenario), scenario := oscen]

ek <- key(eff.dt)
eff.dt <- rbind(eff.dt, refeff.dt) |> setkeyv(ek)

intscns <- eff.dt[, unique(scenario)]
intscn.dt <- scn.dt[scenario %in% intscns][inf_con == FALSE]

intscn.dt[act_alloc == "none", act_alloc := NA]
intscn.dt[pas_alloc == "none", pas_alloc := NA]
intscn.dt[, alloc := fcoalesce(pas_alloc, act_alloc)]

reshaped_wq <- eff.dt[
  intscn.dt, on = .(scenario)
][quar == TRUE][,
  dcast(.SD, outcome + realization + date ~ act_vac, value.var = "cvalue"),
  by = .(alloc, season)
] |> setkey(alloc, season, outcome, realization, date)

reshaped_nq <- eff.dt[
  intscn.dt, on = .(scenario)
][quar == FALSE][,
  dcast(.SD, outcome + realization + date ~ act_vac, value.var = "cvalue"),
  by = .(alloc, season)
] |> setkey(alloc, season, outcome, realization, date)

agg_rank <- function(dt) dt[, c("mscn", "nba", "nbadel") := {
  mat <- .SD |> as.matrix()
  rnk <- (mat |> apply(1, order))
  del <- mat[cbind(1:.N, rnk[2,])] - mat[cbind(1:.N, rnk[1,])]
  .(rnk[1,], rnk[2,], del)
}, .SDcols = !patterns("alloc|season|outcome|realization|date|mscn|nba|nbadel")]

agg_rank(reshaped_wq)
agg_rank(reshaped_nq)

plot_form <- function(dt) {
  res_dt <- dt[, .(
    alloc, season, outcome, realization, date, mscn, nbadel
  )][, refv := nbadel/max(nbadel), by = .(outcome, alloc)]

  reord <- res_dt[
    outcome == "deaths",
    .(refv = sum(refv)),
    by = .(realization, mscn, alloc, season)
  ][,
    .(refv = max(refv), mscn = mscn[which.max(refv)]),
    by = .(alloc, season, realization)
  ][order(mscn, refv), .(realization, neword = 1:.N), by = .(alloc, season)]

  scn_names <- setdiff(names(dt), c(key(dt), "mscn", "nba", "nbadel"))

  res_dt[reord, on = .(realization, alloc, season), neword := neword]
  res_dt[, scenario := factor(scn_names[mscn], levels = scn_names, ordered = TRUE)]

  res_dt
}

wq_plot_dt <- plot_form(reshaped_wq)
nq_plot_dt <- plot_form(reshaped_nq)

saveRDS(wq_plot_dt, gsub("\\.rds$", "_wq.rds", tail(.args, 1)))
saveRDS(nq_plot_dt, gsub("\\.rds$", "_nq.rds", tail(.args, 1)))
