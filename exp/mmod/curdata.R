
# working / invocation directory should be exp/mmod

.args <- if (interactive()) c(
  "smhdata.rds",
  "../active-vac/fig/process/validation.rds",
  "curdata.rds"
) else commandArgs(trailingOnly = TRUE)

library(data.table)

# the submitted SMH data
smh_dt <- readRDS(.args[1])

# the model outputs used as demo time series
# deaths already converted to weekly incidence, but need to convert
# rcase and rhosp to weekly
update_dt <- readRDS(.args[2])[,
  .SD, .SDcols = c("season", "realization", "date", "rcase", "rdeath", "rhosp")
][season != 0]

update_dt$season <- NULL

update_dt[
  order(realization, date),
  c("ccase", "chosp") := .(cumsum(rcase), cumsum(rhosp)), by = realization
]

rev_dt <- (update_dt[
  !is.na(rdeath)
][, .(
  date, case = c(0, diff(ccase)),
  hosp = c(0, diff(chosp)), death = rdeath
), by = realization] |> melt.data.table(
  id.vars = c("realization", "date"), variable.name = "target"
))[between(date, smh_dt[, min(date)], smh_dt[, max(date)])]


saveRDS(rev_dt, tail(.args, 1))