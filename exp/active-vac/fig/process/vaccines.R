.pkgs <- c("data.table", "MMWRweek")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv",
  file.path("fig", "input", "vaccines.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

cdc = fread(.args[2])[
  `Age group` == "all_ages_adj" &
  `Vaccine product` == "all_types" & outcome == "case"
][, date := {
  spl <- lapply(tstrsplit(`MMWR week`,"(?<=.{4})", perl = TRUE), as.integer)
  MMWRweek2Date(spl[[1]], spl[[2]])
}][,
   brkthruRatio := `Vaccinated with outcome`/(`Vaccinated with outcome` + `Unvaccinated with outcome`)
][,
  vaxOutcomeP10k := `Vaccinated with outcome` * (1e4/`Fully vaccinated population`)
][date <= endday]

saveRDS(cdc, tail(.args, 1))
