
# working / invocation directory should be exp/mmod

.args <- if (interactive()) c(
  "smhdata.rds",
  "plot_log21.csv",
  "curdata.rds"
) else commandArgs(trailingOnly = TRUE)

library(data.table)

# the submitted SMH data
smh_dt <- readRDS(.args[1])

# the model outputs used as demo time series
update_dt <- fread(
  .args[2],
  select = c("date", "rcase", "rdeath", "rhosp")
) |> setnames(\(n) gsub("^r", "", n)) |>
  melt.data.table(id.vars = "date", variable.name = "target")

update_dt[,
  target := gsub(".*(case|death|hosp)", "\\1", target) |>
  factor(levels = c("case", "hosp", "death"), ordered = TRUE)
]

dater <- smh_dt[, range(date)]
dater[1] <- dater[1] - 7*4
dater[2] <- dater[2] + 7*4
tardates <- seq(from = dater[1], to = dater[2], by = 7)
update_dt[, cvalue := cumsum(value)]
update_dt <- update_dt[date %in% tardates]
update_dt[, dvalue := c(0, diff(cvalue)), by = target]

update_dt <- update_dt[date != tardates[1]]
update_dt$value <- update_dt$cvalue <- NULL
setnames(update_dt, "dvalue", "value")

saveRDS(update_dt, tail(.args, 1))