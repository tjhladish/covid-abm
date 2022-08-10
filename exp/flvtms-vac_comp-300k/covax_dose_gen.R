#' SOURCES
#' WHO SAGE Prioritized Infectious Disease and Economic Modelling Questions
#' COVID-19 vaccination in Sindh Province, Pakistan: a modelling study of health-impact and cost-effectiveness (SI)

library(data.table)
library(lubridate)

#if (interactive()) { setwd("~/documents/work/covid-abm/exp/flvtms-vac_comp-300k/") }

.args <- if (interactive()) c(
  "active_vax_counterfactual_doses_50k.txt"
) else commandArgs(trailingOnly = TRUE)

dose_file <- fread(.args[1])
dose_file[, n_doses_p10k := 0]

global_pop                   = 7.8e9
covax_participation_fraction = ...
est_covax_pop                = global_pop * covax_participation_fraction

qtr_dates = list(seq.Date(from = ymd("2021-04-01"), to = ymd("2021-06-30"), by = "days"),
                 seq.Date(from = ymd("2021-07-01"), to = ymd("2021-09-30"), by = "days"),
                 seq.Date(from = ymd("2021-10-01"), to = ymd("2021-12-31"), by = "days"),
                 seq.Date(from = ymd("2022-01-01"), to = ymd("2022-03-31"), by = "days"))

qtr_total_doses = c(1e8, 4e8, 6e8, 8e8)

qtr_doses_p10k = qtr_total_doses * 1e4 / est_covax_pop

for (i in 1:length(qtr_dates)) {
    dose_file[is_urg == 1 & bin_min == 5 & dose == 1 & date %in% qtr_dates[[i]], n_doses_p10k := qtr_doses_p10k[i]/length(qtr_dates[[i]])]
}

fwrite(dose_file[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)], file = "covax_doses.txt", sep = " ")
