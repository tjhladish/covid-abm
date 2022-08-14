.pkgs <- c("data.table")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "reporting_dump.csv",
  file.path("fig", "input", "detection.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

.args[2] |> fread() |> saveRDS(file = tail(.args, 1))
