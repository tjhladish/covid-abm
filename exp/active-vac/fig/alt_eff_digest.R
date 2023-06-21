
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "covid-active-v5.1.sqlite"
  )),
  file.path("fig", "process", "alt_eff.rds")
) else commandArgs(trailingOnly = TRUE)

magicdate <- as.Date("2020-12-01")

abcreader <- function(
  pth,
  pmcols = c("serial", "realization", "quar", "pas_vac", "act_vac", "pas_alloc", "act_alloc", "inf_con"),
  metacols = c("serial", "date", "inf", "sev", "deaths"),
  datelim = magicdate,
  reallimit = if (interactive()) 10,
  verbose = interactive()
) {

  wherelim <- if (!is.null(reallimit)) sprintf(" WHERE realization < %i", reallimit) else ""
  stmt <- sprintf(
    "SELECT %s FROM par JOIN met USING(serial)%s;",
    paste(pmcols, collapse=", "),
    wherelim
  )
  if (verbose) warning("Issuing: ", stmt)

  db <- dbConnect(SQLite(), pth)
  dt <- db |> dbGetQuery(stmt) |> as.data.table()
  dt[, realization := as.integer(realization) ] |> setkey(serial)

  if (wherelim != "") {
    srls <- dt[, unique(serial)]
    wherelim <- sprintf(" WHERE serial IN (%s)", paste(srls, collapse = ","))
  }

  mstmt <- sprintf(
    "SELECT %s FROM meta%s;",
    metacols |> paste(collapse=", "),
    wherelim
  )

  if (verbose) warning("Issuing: ", mstmt)
  meta.dt <- db |> dbGetQuery(mstmt) |> as.data.table()

  meta.dt <- meta.dt[, date := as.Date(date) ] |> setkey(serial, date)
  if (!is.null(datelim)) meta.dt <- meta.dt[ date >= datelim ]

  ret <- list(pars = dt, meta = meta.dt)

  if (verbose) { # also fetch out met table
    metstmt <- sprintf(
      "SELECT * FROM met%s;",
      wherelim
    )
    mets.dt <- db |> dbGetQuery(metstmt) |> as.data.table()
    ret$ms <- mets.dt
  }

  dbDisconnect(db)

  return(ret)
}

reserialize <- function(dt, offset = 0L) {
  lapply(dt, function(sdt) sdt[, serial := offset + .GRP, by=serial])
  dt
}

scenarize <- function(dt, offset = 0L) {
  dt$pars[, scenario := offset + 1L:.N, by=.(realization)]
  setkey(dts$pars, serial, realization)
  lapply(dt[-(names(dt) == "pars")], function(sdt) sdt[
    dt$pars, on =.(serial),
    c("scenario", "realization") := .(scenario, realization)
  ])
  dt
}

dts <- head(.args, -1) |> abcreader(datelim = NULL)

reserialize(dts)
scenarize(dts)

genordfac <- \(lvl) return(\(x) lvl[x+1] |> factor(levels = lvl, ordered = TRUE))

funs <- list(
  quar = as.logical,
  pas_vac = as.logical,
  act_vac = genordfac(c("none", "ring", "risk", "age")),
  pas_alloc = genordfac(c("none", "LS", "MS", "HS", "USA")),
  act_alloc = genordfac(c("none", "LS", "MS", "HS", "USA")),
  inf_con = \(x) x == 2,
  scenario = \(x) x
)

scn.dt <- dts$pars[
  realization == 0, .SD, .SDcols = -c("serial", "realization")
][, sapply(names(funs), \(nm) funs[[nm]](.SD[[nm]]), USE.NAMES = TRUE, simplify = FALSE)]
dts$pars <- NULL

#' @examples
#' cmp.dt <- rbind(
#'   dts$meta[, .(cdeath = sum(deaths)*37.5, source = "meta"), by=.(scenario, realization)],
#'   dts$ms[, .(cdeath = tot_cumul_deaths, source = "mets"), by=.(scenario, realization)]
#' )[
#'   ,.(value = median(cdeath)), by=.(scenario, source)
#' ][scn.dt, on=.(scenario)][inf_con == FALSE][,
#'   .(scenario, source, quar, pas_vac, act_vac, alloc = fifelse(pas_vac, pas_alloc, act_alloc), value)
#' ]
#'
#' require(ggplot2); require(ggh4x)
#' ggplot(cmp.dt) + aes(x=alloc, y=value, shape = source, color = act_vac) +
#'   facet_nested("Quar" + quar ~ .) +
#'   geom_point() + theme_minimal() + scale_y_continuous("Tot. Deaths")

meta.dt <- dts$meta[, .SD, .SDcols = -c("serial")] |>
  setkey(scenario, realization, date)

dts$meta <- NULL
gc()

meta.dt <- meta.dt |>
  (\(dt) melt.data.table(dt, key(dt), variable.name = "outcome"))() |>
  setkey(scenario, realization, outcome, date)

meta.dt[, c.value := cumsum(value), by=setdiff(key(meta.dt), "date")]

rm(dts)
gc()

#' the baseline scenarios are:
#'  - non "active" distribution scenarios
#'  - with some passive allocation
#'  - with no additional NPIs (i.e. quarantine program)
basescnid <- scn.dt[
  (act_vac == "none") & (pas_vac == TRUE) & (quar == FALSE)
][,
  unique(scenario)
]

#' for this comparison, also remove the do-nothing scenarios
excludescns <- scn.dt[
  !(scenario %in% basescnid)
][
  ((pas_vac == FALSE) & (act_vac == "none"))
][,
  unique(scenario)
]

ref.dt <- meta.dt[
  (scenario %in% basescnid)# & (outcome != "doses")
][
  scn.dt,
  c("alloc", "inf_con") := .(pas_alloc, inf_con),
  on=.(scenario)
]

int.dt <- meta.dt[
  !(scenario %in% c(basescnid, excludescns))# & (outcome != "doses")
][scn.dt,
  c("alloc", "inf_con", "quar") := .(fifelse(pas_vac, pas_alloc, act_alloc), inf_con, quar),
  on=.(scenario)
]

rm(meta.dt)
gc()

int.dt[
  ref.dt[, .SD, .SDcols = -c("scenario")],
  c("averted", "c.effectiveness") := .(i.value-value, fifelse(c.value == i.c.value, 0, (i.c.value-c.value)/i.c.value)),
  on=.(alloc, inf_con, realization, outcome, date)
]

rm(ref.dt)
gc()

saveRDS(
  int.dt[,c(key(int.dt), "value", "averted", "c.effectiveness"), with = FALSE],
  tail(.args, 1)
)

dts <- head(.args, -1) |> abcreader(metacols = c("serial","date","Rt"))

reserialize(dts)
scenarize(dts)

meta.dt <- dts$meta[, .SD, .SDcols = -c("serial")] |>
  setkey(scenario, realization, date)

saveRDS(meta.dt, gsub("\\.rds","-rt.rds", tail(.args, 1)))

rm(meta.dt)
gc()
