
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "covid-active-v5.0.sqlite"
  )),
  file.path("fig", "process", "digest.rds")
) else commandArgs(trailingOnly = TRUE)

magicdate <- as.Date("2020-12-01")

abcreader <- function(
    pth,
    pmcols = c("serial", "realization", "quar", "pas_vac", "act_vac", "pas_alloc", "act_alloc", "inf_con"),
    metacols = c("serial", "date", "inf", "sev", "deaths", "std_doses + urg_doses AS doses"),
    datelim = magicdate,
    reallimit = if (interactive()) 10,
    verbose = interactive()
) {

  wherelim <- if (!is.null(reallimit)) sprintf(" WHERE realization < %i", reallimit) else ""
  stmt <- sprintf(
    "SELECT %s FROM par JOIN met USING(serial)%s ORDER BY serial;",
    paste(pmcols, collapse=", "),
    wherelim
  )
  if (verbose) warning("Issuing: ", stmt)

  db <- dbConnect(SQLite(), pth)
  dt <- db |> dbGetQuery(stmt) |> as.data.table()
  dt[, realization := as.integer(realization) ]

  if (wherelim != "") {
    srls <- dt[, unique(serial)]
    wherelim <- sprintf(" WHERE serial IN (%s)", paste(srls, collapse = ","))
  }

  mstmt <- sprintf(
    "SELECT %s FROM meta%s ORDER BY serial, date;",
    metacols |> paste(collapse=", "),
    wherelim
  )

  if (verbose) warning("Issuing: ", mstmt)
  meta.dt <- db |> dbGetQuery(mstmt) |> as.data.table()
  meta.dt[, date := as.Date(date) ][ date >= datelim ]

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

dts <- head(.args, -1) |> abcreader()

setkey(dts$pars, serial, realization)
reserialize(dts)
scenarize(dts)

meta.dt <- dts$meta[, .SD, .SDcols = -c("serial")] |>
  setkey(scenario, realization, date)

dts$meta <- NULL
gc()

meta.dt <- meta.dt |>
  (\(dt) melt.data.table(dt, key(dt), variable.name = "outcome"))() |>
  setkey(scenario, realization, outcome, date)

meta.dt[, c.value := cumsum(value), by=setdiff(key(meta.dt), "date")]

scn.dt <- dts$pars[
  realization == 0, .SD, .SDcols = -c("serial", "realization")
]
dts$pars <- NULL
rm(dts)
gc()

# WARNING: MAGIC NUMBER
basescnid <- 1

ref.dt <- meta.dt[(scenario == basescnid) & (outcome != "doses")]
int.dt <- meta.dt[
  (scenario != basescnid) & (outcome != "doses")
]

doses.dt <- meta.dt[outcome == "doses"]
saveRDS(doses.dt, gsub("\\.rds","-doses.rds", tail(.args, 1)))

rm(meta.dt, doses.dt)
gc()

int.dt[
  ref.dt[, .SD, .SDcols = -c("scenario")],
  c("averted", "c.effectiveness") := .(i.value-value, (i.c.value-c.value)/i.c.value),
  on=.(realization, outcome, date)
]

saveRDS(ref.dt, gsub("\\.rds","-ref.rds", tail(.args, 1)))

rm(ref.dt)
gc()

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

saveRDS(
  int.dt[,c(key(int.dt), "value", "averted", "c.effectiveness"), with = FALSE],
  tail(.args, 1)
)

rm(int.dt)
gc()

genordfac <- \(lvl) return(\(x) lvl[x+1] |> factor(levels = lvl, ordered = TRUE))

funs <- list(
  quar = as.logical,
  pas_vac = as.logical,
  act_vac = genordfac(c("none", "ring", "risk")),
  pas_alloc = genordfac(c("none", "LIC", "MIC", "HIC", "USA")),
  act_alloc = genordfac(c("none", "LIC", "MIC", "HIC", "USA")),
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
