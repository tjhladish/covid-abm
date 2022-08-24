.pkgs <- c("data.table")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "population-pseudo-300K.csv",
  "state_based_counterfactual_doses.csv",
  file.path("fig", "input", "seroprev.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

bin_pops = fread(.args[2], select = "age")$age |>
  cut(breaks = c(0, 5, 12, 18, 65, 120), right = F) |>
  table() |> data.table(bin_min = c(0, 5, 12, 18, 65))
poptot <- bin_pops[, sum(N)]

vax = fread(.args[3])[
  (ref_location == "FL") & (is_urg == 0) & (date <= endday)
][
  bin_pops[bin_min != 0, .(bin_min, pop = N)], on=.(bin_min)
][,
  doses := n_doses_p10k * (pop / 1e4)
][,
  .(doses = sum(doses)), keyby = .(date, dose)
][, cov := cumsum(doses)/poptot, by=.(dose) ]

saveRDS(vax, tail(.args, 1))
