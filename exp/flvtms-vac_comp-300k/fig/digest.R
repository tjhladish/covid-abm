
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "covid-active-vac_w-meta.sqlite",
    "covid-active-vac_ring-vac-alternatives_w-meta.sqlite",
    "covid-active-only.sqlite"
  )),
  file.path("fig", "digest.rds")
) else commandArgs(trailingOnly = TRUE)

abcreader <- function(
  pth,
  pmcols = c("serial", "realization", "vac", "dose_file", "active_vax_strat", "quarantine"),
  metacols = c("serial", "date", "inf", "sev", "deaths", "std_doses + urg_doses AS doses")
) {

  db <- dbConnect(SQLite(), pth)
  stmt <- pmcols |> paste(collapse=", ") |> sprintf(
    fmt = "SELECT %s FROM par JOIN met USING(serial) ORDER BY serial;"
  )
  dt <- db |> dbGetQuery(stmt) |> as.data.table()
  dt[, realization := as.integer(realization) ]

  mstmt <- metacols |> paste(collapse=", ") |> sprintf(
    fmt = "SELECT %s FROM meta ORDER BY serial, date;"
  )
  meta.dt <- db |> dbGetQuery(mstmt) |> as.data.table()
  meta.dt[, date := as.Date(date) ]

  dbDisconnect(db)
  return(list(pars = dt, meta = meta.dt))
}

reserialize <- function(dt, offset = 0L) {
  dt$pars[, serial := offset + .GRP, by=serial]
  dt$meta[, serial := offset + .GRP, by=serial]
  dt
}

scenarize <- function(dt, offset = 0L) {
  dt$pars[, scenario := offset + 1L:.N, by=.(realization)]
  dt$meta[
    dt$pars, on =.(serial),
    c("scenario", "realization") := .(scenario, realization)
  ]
  dt
}

metakeys <- c("scenario", "realization", "outcome", "date")

merge.dt <- head(.args, -1) |> lapply(abcreader) |> (\(alldts) Reduce(
  f = \(ldts, rdts) {
    serialoffset <- ldts[[1]]$pars[, max(serial)]
    scenaroffset <- ldts[[1]]$pars[, max(scenario)]
    return(
      rdts |>
      reserialize(offset = serialoffset) |>
      scenarize(offset = scenaroffset) |> list()
    )
  }, x = alldts[-1], init = alldts[[1]] |> reserialize() |> scenarize() |> list(),
  accumulate = TRUE
))() |> (\(alldt) {
  meta <- alldt |>
    lapply(\(pair) pair$meta[
      pair$pars, on=.(serial), c("scenario", "realization") := .(scenario, realization)
    ][, .SD, .SDcols = -c("serial")]) |>
    rbindlist() |>
    melt(
      id.vars = setdiff(metakeys, "outcome"), variable.name = "outcome"
    ) |>
    setkeyv(cols = metakeys)
  meta[, c.value := cumsum(value), by=c(setdiff(metakeys, "date"))]

  scn <- alldt |> lapply(
    \(pair) pair$pars[realization == 0, .SD, .SDcols = -c("serial", "realization")]
  ) |> rbindlist()

  return(list(meta = meta, scn = scn))
})()

caststockpile <- function(x) factor(
  x, levels = c("low", "low-matched", "high", "high-matched"), ordered = TRUE
)

castaction <- function(x) factor(
  x, levels = c("none", "qonly", "vonly", "vandq"), ordered = TRUE
)

castactive <- function(x) factor(
  x, levels = c("none", "passive", "risk", "ring"), ordered = TRUE
)

#' translate serial into meaningful values
scn.dt.base <- dt.base[[1]][realization == 0, .(
  scenario,
  stockpile = caststockpile(c("low", "high")[dose_file]),
  action = castaction(fifelse(
    vac & quarantine, "vandq",
    fifelse(vac == 1, "vonly",
    fifelse(quarantine == 1, "qonly",
      "none"
  )))),
  active = castactive(fifelse(
    active_vax_strat > 0, "ring",
    fifelse(vac == 1, "passive",
    "none"
  )))
)]

scn.dt.alt <- dt.alt[[1]][realization == 0, .(
  scenario,
  stockpile = caststockpile(c("low-matched")),
  action = castaction(fifelse(
    vac & quarantine, "vandq",
    fifelse(vac == 1, "vonly",
    fifelse(quarantine == 1, "qonly",
      "none"
  )))),
  active = castactive(
    fifelse(active_vax_strat > 0, "risk",
    fifelse(vac == 1, "passive",
      "none"
  )))
)]

# no intervention
scn.dt.base[scenario == 1, stockpile := NA ]
# no vaccination, but contact tracing => quarantine
scn.dt.base[scenario == 5, stockpile := NA ]

scn.dt <- rbind(
  scn.dt.base,
  scn.dt.alt
)


gc()

longmeta.dt <- setkey(
  melt(meta.dt, id.vars = c("scenario", "realization", "date"), variable.name = "outcome"),
  scenario, realization, outcome, date
)
#' n.b.: this works because key'ed (ordered) by date
longmeta.dt[,
  c.value := cumsum(value),
  by=.(scenario, realization, outcome)
]

#' no longer need meta.dt
rm(meta.dt)
gc()

#' non-intervention; remove scenario column, so it doesn't appear in joins
#' n.b. this will need to change if there are multiple non-interventions
ref.dt <- longmeta.dt[(scenario == 1) & (outcome != "doses")]
doses.dt <- longmeta.dt[outcome == "doses"]

saveRDS(scn.dt, gsub("\\.rds","-key.rds", tail(.args, 1)))
saveRDS(ref.dt, gsub("\\.rds","-ref.rds", tail(.args, 1)))
saveRDS(doses.dt, gsub("\\.rds","-doses.rds", tail(.args, 1)))

rm(scn.dt, doses.dt)
gc()

#' interventions
int.dt <- setkeyv(
  longmeta.dt[(scenario != 1) & (outcome != "doses")][ref.dt[, .SD, .SDcols = -c("scenario")], on=.(realization, date, outcome)],
  key(longmeta.dt)
)
int.dt[,
  c("averted", "c.averted") := .(i.value-value, i.c.value-c.value)
]

rm(longmeta.dt)
gc()

int.dt[, c.effectiveness := c.averted/i.c.value ]

saveRDS(int.dt, tail(.args, 1))

