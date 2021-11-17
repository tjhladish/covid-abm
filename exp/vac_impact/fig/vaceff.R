
require(data.table)

.args <- if (interactive()) c(
  "rawresults.rds",
  "vaceff.rds"
) else commandArgs(trailingOnly = TRUE)

src <- readRDS(.args[1])[var != "Rt"]

#' contemplates multiple vac states, but assumes 0 == no intervention reference
novac <- src[vac == 0][, -c("serial", "seed", "variable", "week")]
#' keep vac label here to allow for multiple vaccination scenarios to match against
withvac <- src[vac != 0][, -c("serial", "seed", "vac", "variable", "week")]

#' create matched series, compute comparisons
res <- withvac[
  novac,
  on = .(realization, var, date, mutation), nomatch = 0
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