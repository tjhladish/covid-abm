
# working / invocation directory should be exp/mmod

.args <- if (interactive()) c(
  "smhdata.rds"
) else commandArgs(trailingOnly = TRUE)

library(data.table)

smh_url <- "https://raw.githubusercontent.com/midas-network/covid19-scenario-modeling-hub/master/data-processed/UF-ABM/2021-12-29-UF-ABM.csv"

# the submitted SMH data for our model
smh_dt <- fread(
  smh_url,
  select = c(
    "scenario_name", "target", "target_end_date", "quantile", "type", "value"
  )
)[
  (type == "quantile") & (quantile %in% c(0.025, 0.5, 0.975))
][
  target %like% "^.* inc .*$"
][, target := gsub(".*(case|death|hosp)", "\\1", target) |>
    factor(levels = c("case", "hosp", "death"), ordered = TRUE)
] |> 
  setnames(c("target_end_date", "scenario_name"), c("date", "scenario")) |>
  setkey(scenario, target, quantile, date) |> setcolorder()

# n.b. per SMH meta data, our model was 375K individuals,
# but our outputs were re-scaled to SMH population
popFL <- 21477737

# for comparison, translate to per 10K
smh_dt[, value := value / (popFL / 10000) ]

saveRDS(smh_dt, tail(.args, 1))
