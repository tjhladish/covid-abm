
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "covid-active-v3.1.sqlite"
  )),
  file.path("fig", "process", "alt_eff.rds")
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

genordfac <- \(lvl) return(\(x) lvl[x+1] |> factor(levels = lvl, ordered = TRUE))

funs <- list(
  quar = as.logical,
  pas_vac = as.logical,
  act_vac = genordfac(c("none", "ring", "risk")),
  pas_alloc = genordfac(c("none", "FL", "FL+ring", "COVAX", "MIC")),
  act_alloc = genordfac(c("none", "ring", "ringmonth", "COVAX", "MIC")),
  inf_con = \(x) x == 2
)

scn.dt <- dts$pars[
  realization == 0, .SD, .SDcols = -c("serial", "realization")
][, sapply(names(funs), \(nm) funs[[nm]](.SD[[nm]]), USE.NAMES = TRUE, simplify = FALSE)][, scenario := 1:.N ]
dts$pars <- NULL
rm(dts)
gc()

# WARNING: MAGIC NUMBER
basescnid <- scn.dt[
  (quar == FALSE) & (act_vac == "none") & (pas_alloc %in% c("COVAX", "MIC", "FL+ring")),
  unique(scenario)
]

nomatchscns <- scn.dt[
  !(scenario %in% basescnid)
][
  ((pas_alloc == "FL" & act_alloc == "none") | (pas_alloc == "none" & act_alloc == "none")),
  unique(scenario)
]

ref.dt <- meta.dt[
  (scenario %in% basescnid) & (outcome != "doses")
][scn.dt, c("alloc", "inf_con") := .(pas_alloc, inf_con), on=.(scenario)]
int.dt <- meta.dt[
  !(scenario %in% c(basescnid, nomatchscns)) & (outcome != "doses")
][scn.dt, c("alloc", "inf_con") := .(
  fifelse(pas_alloc == "none",
    (as.integer(act_alloc) - 1L) |> funs$pas_alloc(),
    fifelse(pas_alloc == "FL", funs$pas_alloc(2), pas_alloc)
  ),
  inf_con
), on=.(scenario)]

rm(meta.dt)
gc()

int.dt[
  ref.dt[, .SD, .SDcols = -c("scenario")],
  c("averted", "c.effectiveness") := .(i.value-value, (i.c.value-c.value)/i.c.value),
  on=.(alloc, inf_con, realization, outcome, date)
]

rm(ref.dt)
gc()

saveRDS(
  int.dt[,c(key(int.dt), "value", "averted", "c.effectiveness"), with = FALSE],
  tail(.args, 1)
)
