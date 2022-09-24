.pkgs <- c("data.table")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "reporting_by_day.txt",
  file.path("fig", "process", "detection.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

.args[2] |> fread() |>
  DT(,
     mild := asymp + mild - asymp*mild
  ) |>
  DT(,
     severe := mild + severe - mild*severe
  ) |>
  DT(,
     crit := severe + crit - severe*crit
  ) |>
  DT(,
     death := crit + death - crit*death
  ) |>
  saveRDS(file = tail(.args, 1))
