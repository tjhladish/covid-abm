.pkgs <- c("data.table", "MMWRweek")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv",
  file.path("fig", "input", "hospitals.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

per10k <- 1e4/pop["florida"]

hhsHosp = fread(.args[2])[order(date)][
  state == 'FL'
][,
  date := as.Date(date)
][,
  hospInc := previous_day_admission_adult_covid_confirmed * per10k
][date <= endday]

saveRDS(hhsHosp, tail(.args, 1))
