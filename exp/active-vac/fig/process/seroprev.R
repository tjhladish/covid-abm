.pkgs <- c("data.table")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "CDC_seroprev_long.csv",
  file.path("fig", "input", "seroprev.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

seroprev = fread(.args[2])[
  !is.na(est)
][,
  unique(.SD), keyby=date
][, run := {
  res <- rle(est)
  rep(1:length(res$lengths), res$lengths)
}][,.(
  start = date[1], end = date[.N],
  lower = lower[1], est = est[1], upper = upper[1]
), by=run
]

saveRDS(seroprev, tail(.args, 1))
