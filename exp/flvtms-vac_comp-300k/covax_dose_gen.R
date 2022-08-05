#' SOURCES
#' WHO SAGE Prioritized Infectious Disease and Economic Modelling Questions
#' COVID-19 vaccination in Sindh Province, Pakistan: a modelling study of health-impact and cost-effectiveness (SI)

library(data.table)
library(lubridate)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/flvtms-vac_comp-300k/") }

.args <- if (interactive()) c(
  "active_vax_counterfactual_doses.txt"
) else commandArgs(trailingOnly = TRUE)

dose_file <- fread(.args[1])
dose_file[, n_doses_p10k := 0]

total_num_countries = 195
covax_participants = 166
covax_participation_fraction = covax_participants/total_num_countries

global_pop <- 7.8e9
est_covax_pop <- global_pop * covax_participation_fraction

q2_2021_dates <- seq.Date(from = ymd("2021-04-01"), to = ymd("2021-06-30"), by = "days")
q3_2021_dates <- seq.Date(from = ymd("2021-07-01"), to = ymd("2021-09-30"), by = "days")
q4_2021_dates <- seq.Date(from = ymd("2021-10-01"), to = ymd("2021-12-31"), by = "days")
q1_2022_dates <- seq.Date(from = ymd("2022-01-01"), to = ymd("2022-03-31"), by = "days")

q2_2021_total_doses <- 1e8
q3_2021_total_doses <- 4e8
q4_2021_total_doses <- 6e8
q1_2022_total_doses <- 8e8

q2_2021_doses_p10k <- q2_2021_total_doses * (1e4/ est_covax_pop)
q3_2021_doses_p10k <- q3_2021_total_doses * (1e4/ est_covax_pop)
q4_2021_doses_p10k <- q4_2021_total_doses * (1e4/ est_covax_pop)
q1_2022_doses_p10k <- q1_2022_total_doses * (1e4/ est_covax_pop)

dose_file[is_urg == 1 & bin_min == 5 & date %in% q2_2021_dates, n_doses_p10k := q2_2021_doses_p10k]
dose_file[is_urg == 1 & bin_min == 5 & date %in% q3_2021_dates, n_doses_p10k := q3_2021_doses_p10k]
dose_file[is_urg == 1 & bin_min == 5 & date %in% q4_2021_dates, n_doses_p10k := q4_2021_doses_p10k]
dose_file[is_urg == 1 & bin_min == 5 & date %in% q1_2022_dates, n_doses_p10k := q1_2022_doses_p10k]

fwrite(dose_file[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)], file = "covax_doses.txt", sep = " ")