#!/usr/bin/env Rscript
library(lubridate)
library(tidyr)
args = commandArgs(trailingOnly=TRUE)

data_dir = '.'
if (length(args) == 1) {
  data_dir = args[1]
} else {
  stop("Pass data directory as command-line argument")
}

### CHANGE THESE PARAMETERS TO CONTROL THE PROGRAM
num_replicates = 100
percentiles = c(0.00, 0.01, 0.025, seq(0.05, 0.95, by=0.05), 0.975, 0.99, 1.00)
scen_ids = c('A-2022-01-09', 'B-2022-01-09', 'C-2022-01-09', 'D-2022-01-09')
scen_names = c('optSev_highIE', 'optSev_lowIE', 'pessSev_highIE', 'pessSev_lowIE')
model_projection_date = '2022-01-13'
location = 12
fl_pop = 21538187

# files named plot_log#.csv where # is the serial
tmp = read.table(file=grep(paste0("^", data_dir, "/plot_log(\\d+).csv"), list.files(data_dir, full.names=T), value=T)[1],
                 sep=',', stringsAsFactors=F, header=T)
d1 = data.frame(epiweek = integer(), inc_rcase = double(), inc_rdeath = double(), inc_hosp = double(),
                cuml_rcase = double(), cuml_rdeath = double(), cuml_hosp = double(), epiweek_end_date = as.Date(character()))
d2 = data.frame(epiweek = integer(), inc_rcase = double(), inc_rdeath = double(), inc_hosp = double(),
                cuml_rcase = double(), cuml_rdeath = double(), cuml_hosp = double(), epiweek_end_date = as.Date(character()))
d3 = data.frame(epiweek = integer(), inc_rcase = double(), inc_rdeath = double(), inc_hosp = double(),
                cuml_rcase = double(), cuml_rdeath = double(), cuml_hosp = double(), epiweek_end_date = as.Date(character()))
d4 = data.frame(epiweek = integer(), inc_rcase = double(), inc_rdeath = double(), inc_hosp = double(),
                cuml_rcase = double(), cuml_rdeath = double(), cuml_hosp = double(), epiweek_end_date = as.Date(character()))

scen2_min = num_replicates * 1
scen3_min = num_replicates * 2
scen4_min = num_replicates * 3

for (file in list.files(data_dir, full.names=T)) {
  if (grepl(paste0("^", data_dir, "/plot_log(\\d+).csv"), file)) {
    tmp = read.table(file=file, sep=',', stringsAsFactors=F, header=T)
    tmp_serial = as.integer(sub(paste0("^", data_dir, "/plot_log(\\d+).csv"), '\\1', file))
    tmp = tmp[(tmp$date >= "2021-12-19" & tmp$date <= "2022-05-12"),]
    
    cut_tmp = data.frame(inc_rcase = tmp$rcase, inc_rdeath = tmp$rdeath, inc_hosp = tmp$hospInc)
    cut_tmp$epiweek_end_date = ceiling_date(ymd(tmp$date), 'week') - 1
    cut_tmp$epiweek = epiweek(tmp$date)
    cut_tmp$cuml_rcase = cumsum(cut_tmp$inc_rcase)
    cut_tmp$cuml_rdeath = cumsum(cut_tmp$inc_rdeath)
    cut_tmp$cuml_hosp = cumsum(cut_tmp$inc_hosp)
    
    agg_tmp = aggregate(cut_tmp[,-c(4,5)], by=list(epiweek=cut_tmp$epiweek), FUN=sum)
    for (ew in agg_tmp$epiweek) {
      agg_tmp$epiweek_end_date[agg_tmp$epiweek == ew] = as.character(unique(cut_tmp$epiweek_end_date[cut_tmp$epiweek == ew]))
    }
    agg_tmp$epiweek_end_date = as.Date(agg_tmp$epiweek_end_date)
    
    if (tmp_serial < scen2_min) {
      d1 = rbind(d1, agg_tmp)
    } else if (tmp_serial >= scen2_min & tmp_serial < scen3_min) {
      d2 = rbind(d2, agg_tmp)
    } else if (tmp_serial >= scen3_min & tmp_serial < scen4_min) {
      d3 = rbind(d3, agg_tmp)
    } else if (tmp_serial >= scen4_min) {
      d4 = rbind(d4, agg_tmp)
    }
  }
}

pctl_cols = paste0('', percentiles)
scen_d = list(d1, d2, d3, d4)

final_df = data.frame(scenario_id = integer(1), met = character(1), epiweek = integer(1), epiweek_end_date = as.Date('2020-01-01'))
final_df[pctl_cols] = NA
final_colnames = colnames(final_df)
final_df = final_df[-1,]
colnames(final_df) = final_colnames

for (df in 1:length(scen_d)) {
  for (col in colnames(scen_d[[df]])[-c(1,8)]) {
    for (epwk in unique(scen_d[[df]]$epiweek)) {
      tmp_quantiles = quantile(scen_d[[df]][scen_d[[df]]$epiweek==epwk,col], percentiles)
      tmp_epwk_end_date = as.Date(unique(scen_d[[df]]$epiweek_end_date[scen_d[[df]]$epiweek == epwk]))
      
      tmp_bind = data.frame(scenario_id = df, met = col, epiweek = epwk, epiweek_end_date = tmp_epwk_end_date)
      tmp_bind[pctl_cols] = tmp_quantiles
      colnames(tmp_bind) = colnames(final_df)
      final_df = rbind(final_df, tmp_bind)
    }
  }
}

final_df = gather(final_df, quantile, value, '0':'1')
final_df$scenario_name = scen_names[final_df$scenario_id]
final_df$scenario_id = scen_ids[final_df$scenario_id]
final_df$model_projection_date = as.Date(model_projection_date)
final_df$location = location
final_df$type = 'quantile'
final_df$value = final_df$value * (fl_pop/1e4)

sorted_end_dates = sort(as.Date(unique(final_df$epiweek_end_date)))
final_df$target = paste0(match(as.Date(final_df$epiweek_end_date), sorted_end_dates), ' wk ahead ', final_df$met)

tmp = final_df
final_df = data.frame(scenario_id = tmp$scenario_id, scenario_name = tmp$scenario_name, model_projection_date = tmp$model_projection_date,
                      target = tmp$target, target_end_date = tmp$epiweek_end_date, quantile = tmp$quantile, type = tmp$type, location = tmp$location,
                      value = tmp$value)

medians = final_df[final_df$quantile == 0.5,]
medians$type = 'point'
medians$quantile = 'NA'
final_df = rbind(final_df, medians)

final_filename = paste0(model_projection_date, '-UF-ABM.csv')
write.table(final_df, file=final_filename, sep=',', quote=F, row.names=F)
