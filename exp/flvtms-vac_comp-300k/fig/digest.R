
.pkgs <- c("data.table", "RSQLite", "jsonlite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path(c(
    "abc-active-vac.json", "covid-active-vac_w-meta.sqlite"
  )),
  file.path("fig", "digest.rds")
) else commandArgs(trailingOnly = TRUE)

js <- read_json(.args[1])

db <- dbConnect(SQLite(), .args[2])
stmt <- paste0("SELECT P.*, M.* FROM par AS P JOIN met AS M USING(serial);")
dt <- as.data.table(dbGetQuery(db, stmt))
meta.dt <- as.data.table(dbGetQuery(db, "SELECT * FROM meta;"))
dbDisconnect(db)

#' assign scenarios; assert: abc inputs stable orderings of
#' the other pseudo parameter settings when viewing table per realization
dt[, scenario := 1:.N, by=.(realization)]

serial.dt <- dt[,.(serial, scenario, realization, vac, mutation, dose_file, active_vax_strat, quarantine)]

meta.dt[, date := as.Date(date) ]
meta.dt[
  serial.dt, on =.(serial),
  c("scenario", "realization") := .(scenario, realization)
]

#' no longer need serial numbers
meta.dt$serial <- NULL
serial.dt$serial <- NULL

longmeta.dt <- setkey(
  melt(meta.dt, id.vars = c("scenario", "realization", "date"), variable.name = "outcome"),
  scenario, realization, outcome, date
)
#' n.b.: this works because key'ed (ordered) by date
longmeta.dt[, c.value := cumsum(value), by=.(scenario, realization, outcome)]

#' non-intervention; remove scenario column, so it doesn't appear in joins
#' n.b. this will need to change if there are multiple non-interventions
ref.dt <- longmeta.dt[scenario == 1]
ref.dt$scenario <- NULL

#' interventions
int.dt <- setkeyv(
  longmeta.dt[ref.dt, on=.(realization, date, outcome)],
  key(longmeta.dt)
)
int.dt[, c("averted", "c.averted") := .(i.value-value, i.c.value-c.value) ]
int.dt[, c.effectiveness := c.averted/i.c.value ]

#' translate serial into meaningful values
scn.dt <- serial.dt[realization == 0, .(
  scenario,
  stockpile = factor(
    c("low", "high")[dose_file], levels = c("low", "high"), ordered = TRUE
  ),
  action = fifelse(vac & quarantine, "vandq",
           fifelse(vac == 1, "vonly",
           fifelse(quarantine == 1, "qonly",
           "none"
  ))),
  active = fifelse(active_vax_strat > 0, "ring",
           fifelse(vac == 1, "passive",
           "none"
  ))
)]

# no intervention
scn.dt[scenario == 1, stockpile := NA ]
# no vaccination, but contact tracing => quarantine
scn.dt[scenario == 5, stockpile := NA ]

saveRDS(int.dt, tail(.args, 1))
saveRDS(scn.dt, gsub("\\.rds","-key.rds", tail(.args, 1)))
