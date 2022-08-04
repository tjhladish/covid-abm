
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("covid-active-only.sqlite"),
  file.path("fig", "digest-key.rds"),
  file.path("fig", "digest-ref.rds"),
  file.path("fig", "digest-ext.rds")
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

dt.only <- abcreader(.args[1])
dt.key <- readRDS(.args[2])

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

reserialize(dt.only)
scenarize(dt.only, offset = dt.key[, max(scenario)])

caststockpile <- function(x) factor(
  x, levels = c("low", "low-matched", "high", "high-matched", "nomass"), ordered = TRUE
)

castaction <- function(x) factor(
  x, levels = c("none", "qonly", "vonly", "vandq"), ordered = TRUE
)

castactive <- function(x) factor(
  x, levels = c("none", "passive", "risk", "ring"), ordered = TRUE
)

#' translate serial into meaningful values
scn.dt <- dt.only[[1]][realization == 0, .(
  scenario,
  stockpile = caststockpile(c("low", "high", "", "", "low")[dose_file]),
  action = castaction(fifelse(
    vac & quarantine, "vandq",
    fifelse(vac == 1, "vonly",
    fifelse(quarantine == 1, "qonly",
      "none"
  )))),
  active = castactive(fifelse(
    active_vax_strat == 2, "ring",
    fifelse(vac == 1, "risk",
    "none"
  )))
)]

serial.dt <- rbind(dt.only[[1]])[,.(
  serial, scenario, realization
)]

meta.dt <- rbind(dt.only[[2]])[
  serial.dt, on =.(serial),
  c("scenario", "realization") := .(scenario, realization)
]

#' no longer need serial numbers
meta.dt$serial <- NULL
rm(serial.dt, dt.only)
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
doses.dt <- longmeta.dt[outcome == "doses"]

saveRDS(doses.dt, gsub("\\.rds","-doses.rds", tail(.args, 1)))
saveRDS(scn.dt, gsub("\\.rds","-key.rds", tail(.args, 1)))

rm(scn.dt, doses.dt)
gc()

ref.dt <- readRDS(.args[3])

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

