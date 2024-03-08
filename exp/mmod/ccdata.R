
# working / invocation directory should be exp/mmod

.args <- if (interactive()) c(
  "smhdata.rds",
  "ccdata.rds"
) else commandArgs(trailingOnly = TRUE)

library(data.table)
library(covidcast)

# the submitted SMH data
smh_dt <- readRDS(.args[1])

dater <- smh_dt[, range(date)]
dater[1] <- dater[1] - 7*4
dater[2] <- dater[2] + 7*4
tardates <- seq(from = dater[1], to = dater[2], by = 7)

# covidcast values in per 100K, not per 10K => so divide value by 10
ref_dt <- rbind(
  (covidcast::covidcast_signal(
    "jhu-csse",
    signal = "confirmed_cumulative_prop", 
    geo_type = "state", geo_values = "fl",
    start_day = tardates[1], end_day = tail(tardates, 1)
  ) |> setDT())[, .(
    date = time_value, target = "case", value = value/10
  )],
  (covidcast::covidcast_signal(
    "hhs",
    signal = "confirmed_admissions_covid_1d_prop", 
    geo_type = "state", geo_values = "fl",
    start_day = tardates[1], end_day = tail(tardates, 1)
  ) |> setDT())[, .(
    date = time_value, target = "hosp", value = cumsum(value/10)
  )],
  (covidcast::covidcast_signal(
    "jhu-csse",
    signal = "deaths_cumulative_prop", 
    geo_type = "state", geo_values = "fl",
    start_day = tardates[1], end_day = tail(tardates, 1)
  ) |> setDT())[, .(
    date = time_value, target = "death", value = value/10
  )]
)

ref_dt[, target := factor(target, levels = c("case", "hosp", "death"), ordered = TRUE)]

ref_dt <- ref_dt[date %in% tardates]
ref_dt[, dvalue := c(0, diff(value)), by = target]

ref_dt <- ref_dt[date != tardates[1]]
ref_dt$value <- NULL
setnames(ref_dt, "dvalue", "value")
ref_dt[, date := as.IDate(date)]

excessmort <- covidcast::covidcast_signals(
  "nchs-mortality", geo_type = "state", geo_values = "fl",
  signal = c("deaths_allcause_incidence_prop", "deaths_percent_of_expected"),
  time_type = "week", start_day = tardates[1], end_day = tail(tardates, 1)
) |> rbindlist() |> dcast.data.table(time_value ~ signal, value.var = "value")

excessmort[, value := (deaths_percent_of_expected - 100)/100*deaths_allcause_incidence_prop/10 ]

ref_dt <- rbind(ref_dt, excessmort[, .(date = as.IDate(time_value), target = "exdeath", value)])


saveRDS(ref_dt, tail(.args, 1))
