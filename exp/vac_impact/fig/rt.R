
require(data.table)

.args <- if (interactive()) c(
  "rawresults.rds",
  "rt.rds"
) else commandArgs(trailingOnly = TRUE)

rt.dt <- readRDS(.args[1])[var == "Rt"][, .(
  realization, date, var = "c", variable = "Rt", value, intervention = as.logical(vac)
)]

saveRDS(rt.dt, tail(.args, 1))