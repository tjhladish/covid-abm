
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "covid-active-vac_w-meta.sqlite",
    "covid-active-vac_ring-vac-alternatives_w-meta.sqlite"
  )),
  file.path("fig", "digest.rds")
) else commandArgs(trailingOnly = TRUE)

abcreader <- function(
  pth,
  pmcols = c("serial", "realization", "vac", "dose_file", "active_vax_strat", "quarantine")
) {
  db <- dbConnect(SQLite(), pth)
  stmt <- sprintf("SELECT %s FROM par AS P JOIN met AS M USING(serial) ORDER BY serial;", paste(pmcols, collapse = ", "))
  dt <- as.data.table(dbGetQuery(db, stmt))
  dt[, realization := as.integer(realization) ]
  meta.dt <- as.data.table(dbGetQuery(db, "SELECT * FROM meta ORDER BY serial, date;"))[, date := as.Date(date) ]
  meta.dt[, doses := std_doses + urg_doses ]
  meta.dt$std_doses <- meta.dt$urg_doses <- NULL
  dbDisconnect(db)
  return(list(dt, meta.dt))
}

dt.base <- abcreader(.args[1])
dt.alt <- abcreader(.args[2])

reserialize <- function(dt, offset = 0L) {
  dt[[1]][, serial := offset + .GRP, by=serial]
  dt[[2]][, serial := offset + .GRP, by=serial]
}

scenarize <- function(dt, offset = 0L) {
  dt[[1]][, scenario := offset + 1L:.N, by=.(realization)]
  dt[[2]][
    dt[[1]], on =.(serial),
    c("scenario", "realization") := .(scenario, realization)
  ]
}

reserialize(dt.base)
reserialize(dt.alt, offset = dt.base[[1]][, max(serial)])
scenarize(dt.base)
scenarize(dt.alt, offset = dt.base[[1]][, max(scenario)])

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

serial.dt <- rbind(dt.base[[1]], dt.alt[[1]])[,.(
  serial, scenario, realization
)]

meta.dt <- rbind(dt.base[[2]], dt.alt[[2]])[
  serial.dt, on =.(serial),
  c("scenario", "realization") := .(scenario, realization)
]

#' no longer need serial numbers
meta.dt$serial <- NULL
rm(serial.dt, dt.alt, dt.base, scn.dt.base, scn.dt.alt)
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

