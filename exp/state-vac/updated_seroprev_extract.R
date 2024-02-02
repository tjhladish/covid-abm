library(data.table)
library(lubridate)

sero_in = fread("cdc_covid-19_seroprev_data.csv")
fl_sero = sero_in[Site == "FL", .(site               = Site,
                                  date_range         = `Date Range of Specimen Collection`,
                                  ped_seroprev_point = `Rate (%) [0-17 Years Prevalence]`,
                                  ped_seroprev_lower = `Lower CI [0-17 Years Prevalence]`,
                                  ped_seroprev_upper = `Upper CI [0-17 Years Prevalence]`,
                                  all_seroprev_point = `Rate (%) [All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                  all_seroprev_lower = `Lower CI [All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                  all_seroprev_upper = `Upper CI [All Ages Cumulative Prevalence, Rounds 1-30 only]`)]

for (row in 1:nrow(fl_sero)) {
  date_str = fl_sero[row, date_range]
  date1 = mdy(paste0(strsplit(gsub("[[:space:]]", "", date_str), split = "-")[[1]][1], ",", strsplit(strsplit(gsub("[[:space:]]", "", date_str), split = "-")[[1]][2], split = ",")[[1]][2]))
  date2 = mdy(paste0(strsplit(strsplit(gsub("[[:space:]]", "", date_str), split = "-")[[1]][2], split = ",")[[1]][1], ",", strsplit(strsplit(gsub("[[:space:]]", "", date_str), split = "-")[[1]][2], split = ",")[[1]][2]))
  fl_sero[row, `:=` (
    date_start = date1,
    date_end   = date2
  )]
}

fl_sero[11, date_start := ymd("2020-12-28")]

sero_long = data.table()
for (row in 1:nrow(fl_sero)) {
  tmp = data.table()
  tmp[, `:=` (
    date = seq.Date(from = fl_sero[row, date_start], to = fl_sero[row, date_end], by = "day"),
    ped_seroprev_point = fl_sero[row, ped_seroprev_point],
    ped_seroprev_lower = fl_sero[row, ped_seroprev_lower],
    ped_seroprev_upper = fl_sero[row, ped_seroprev_upper],
    all_seroprev_point = fl_sero[row, all_seroprev_point],
    all_seroprev_lower = fl_sero[row, all_seroprev_lower],
    all_seroprev_upper = fl_sero[row, all_seroprev_upper]
  )]
  
  sero_long = rbindlist(list(sero_long, tmp))
}

fwrite(x = sero_long, file = "fl_seroprev_cdc.csv")
