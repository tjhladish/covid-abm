
require(data.table)

.args <- if (interactive()) c(
  "rawresults.rds",
  "vareff.rds"
) else commandArgs(trailingOnly = TRUE)

src <- readRDS(.args[1])[var != "Rt"]

nomut <- src[mutation == 0][, -c("serial", "seed", "mutation", "variable", "week")]
withmut <- src[mutation != 0][, -c("serial", "seed", "variable", "week")]

res <- withmut[
  nomut,
  on = .(realization, var, date, vac), nomatch = 0
]

res[order(date), c(
  "cvalue", "i.cvalue", "averted", "caverted", "eff", "ceff"
) := {
  cv <- cumsum(value)
  icv <- cumsum(i.value)
  av <- i.value - value
  .(cv, icv, av, icv - cv, av/i.value, (icv-cv)/icv)
}, by=.(realization, var, mutation, vac)]

setkey(res, vac, mutation, realization, var, date)

saveRDS(res, tail(.args, 1))