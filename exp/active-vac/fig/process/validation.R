
.pkgs <- c("data.table", "RSQLite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "plotlog"),
  file.path("~", "Downloads", "validation_runs.sqlite"),
  file.path("fig", "process", "validation.rds")
) else commandArgs(trailingOnly = TRUE)

fls <- list.files(.args[2], full.names = TRUE)

extractor <- function(fl) {
  ret <- setnames(fread(fl)[,
    c("crcase", "crdeath") := .(cumsum(rcase), cumsum(rdeath))
  ][,
    VESAvg := filter(VES, rep(1/7, 7), sides = 2)
  ][,
    brkthruRatioAvg := filter(brkthruRatio, rep(1/7, 7), sides = 2)
  ][,
    c("tot_std_doses", "tot_urg_doses") := .(cumsum(std_doses), cumsum(urg_doses))
  ], "serial", "realization")

  # aggregate deaths to week => empirical data is by week
  ret[wday(date) != 7, crdeath := NA]
  ret[, rdeath := NA ]
  ret[wday(date) == 7, rdeath := c(crdeath[1], diff(crdeath)) ]

  # breakthru ratio also on weekly basis, though apparently different wday?
  ret[, brkthru := {
    rel <- (brkthruRatio * rcase)
    rel[is.na(rel)] <- 0
    rel <- cumsum(rel)[wday(date) == 1]
    # breakthrough incidence
    rel <- c(rel[1], diff(rel))
    allinc <- crcase[wday(date) == 1]
    allinc <- c(allinc[1], diff(allinc))
    tmp <- rep(NA, .N)
    tmp[wday(date) == 1] <- rel/allinc
    tmp
  } ]

  ret
}

dt <- rbindlist(lapply(fls, extractor))

db <- dbConnect(SQLite(), .args[3])
seas_dt <- setDT(dbGetQuery(db, "SELECT serial, realization, seed, season FROM par;"))
seas_dt[, seed := bit64::as.integer64(seed) ]
dbDisconnect(db)

dt[seas_dt, on = .(realization = seed), c("realization", "season") := .(i.realization, season)]
dt[, realization := as.integer(realization)]

saveRDS(dt, tail(.args, 1))

