
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "covid-active-v3.0.sqlite"
  )),
  file.path("fig", "process", "digest.rds")
) else commandArgs(trailingOnly = TRUE)

abcreader <- function(
  pth,
  pmcols = c("serial", "realization", "quar", "pas_vac", "act_vac", "pas_alloc", "act_alloc", "inf_con", "ppb_ctrl"),
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

# WARNING: MAGIC NUMBER
basescnid <- 13

ref.dt <- merge.dt$meta[(scenario == basescnid) & (outcome != "doses")]
int.dt <- merge.dt$meta[
  (scenario != basescnid) & (outcome != "doses")
][
  ref.dt[, .SD, .SDcols = -c("scenario")], on=.(realization, date, outcome)
]

scn.dt <- merge.dt$scn
doses.dt <- merge.dt$meta[outcome == "doses"]

rm(merge.dt)

saveRDS(ref.dt, gsub("\\.rds","-ref.rds", tail(.args, 1)))
saveRDS(doses.dt, gsub("\\.rds","-doses.rds", tail(.args, 1)))

rm(ref.dt, doses.dt)

gc()

int.dt[,
  c("averted", "c.averted") := .(i.value-value, i.c.value-c.value)
][,
  c.effectiveness := c.averted/i.c.value
]

#' in general, plotting value, averted, & c.eff
#' potential also want c0.eff = c.eff, but
#' calculated from a later baseline. to compute that, need
#' need c.averted + c.i.value, from new baseline
#' c.i.value = cumsum(value + averted)
#' so, data to store:
#'  * keying values (scenario, realization, date, outcome)
#'  * values (value, averted, [can easily reconstitute c.value, c.averted])
#'  * c.effectiveness
#'  ... then can as necessary reconstitute:
#'  c.value = cumsum(value), c.averted = cumsum(averted)
#'  c0.effectiveness = c.averted (rebased) / cumsum(value + averted) (also rebased)

saveRDS(int.dt[,.(value, averted, c.effectiveness), keyby=.(scenario, realization, outcome, date)], tail(.args, 1))

rm(int.dt)
gc()

genordfac <- \(lvl) return(\(x) lvl[x+1] |> factor(levels = lvl, ordered = TRUE))

funs <- list(
  quar = as.logical,
  pas_vac = as.logical,
  act_vac = genordfac(c("none", "ring", "risk")),
  pas_alloc = genordfac(c("none", "FL", "FL+ring", "COVAX", "MIC")),
  act_alloc = genordfac(c("none", "ring", "ringmonth", "COVAX", "MIC")),
  inf_con = \(x) x == 2
)

res <- scn.dt[, sapply(names(funs), \(nm) funs[[nm]](.SD[[nm]]), USE.NAMES = TRUE, simplify = FALSE)][, scenario := 1:.N ]

saveRDS(res, gsub("\\.rds","-key.rds", tail(.args, 1)))

# WARNING: SERIOUS MAGIC HERE
# TODO refactor when updated scenario runs
# scn.dt[,
#        c("stockpile", "action", "active") := .(
#          caststockpile(dose_file, vac),
#          castaction(vac, quarantine),
#          castactive(active_vax_strat, vac)
#        )
# ]


# caststockpile <- function(df, v) factor(
#   fifelse(v == 0, "none", fifelse(
#     df == 1, "low", fifelse(
#       df == 1.5, "low-matched", fifelse(
#         df == 2, "high", fifelse(
#           df == 5, "covax", "ERROR"
#         ))))), levels = c(
#           "none", "low", "low-matched", "high", "covax"
#         ), ordered = TRUE
# )
#
# castaction <- function(v,q) factor(
#   fifelse(
#     v == q, fifelse(v==1,"vandq", "none"), fifelse(v==1, "vonly", "qonly")
#   ), levels = c("none", "qonly", "vonly", "vandq"), ordered = TRUE
# )
#
# castactive <- function(strat, v) factor(
#   fifelse(
#     strat == 0, fifelse(v==0, "none", "passive"), fifelse(
#       strat == 2, "ring", fifelse(
#         strat == 6, "risk", "ERROR"
#       ))), levels = c("none", "passive", "risk", "ring"), ordered = TRUE
# )
