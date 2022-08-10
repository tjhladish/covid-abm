library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)
library(DBI)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/fl_vac_counterfactual/") }

.args <- if (interactive()) c(
  "ACS_2019_pop_data.csv",
  "cdc_covid-19_vax_data.csv",
  "./dose_data",
  "sim_data_0.sqlite"
) else commandArgs(trailingOnly = TRUE)

db_path <- .args[4]

serial <- sub("sim_data_([[:digit:]]*).sqlite", "\\1", db_path)
figPath <- paste0("./fig/sim_", serial, "_diagnostics")
dir.create(figPath, showWarnings = FALSE)

source("./vax_adj_v7.R")

counterfact_scen_to_ref_loc <- function(scen) {
  switch(scen,
         '0' = { ref_loc = "FL" },
         '1' = { ref_loc = "VT" },
         '2' = { ref_loc = "MS" })
  return(ref_loc)
}

db <- dbConnect(RSQLite::SQLite(), "./covid_counterfactuals_v1.0.sqlite")
stmt <- paste0("select counterfact_scen from par where serial = ", serial, ";")
model_ref_loc <- counterfact_scen_to_ref_loc(as.character(dbGetQuery(db, stmt)))
dbDisconnect(db)


db <- dbConnect(RSQLite::SQLite(), db_path)
stmt <- "select * from vaccination_history;"
vax_hist <- dbGetQuery(db, stmt)
vax_hist <- setDT(vax_hist)
vax_hist[, vax_date := ymd(vax_date)]

stmt <- "select * from age_bins;"
bin_pops <- dbGetQuery(db, stmt)
bin_pops <- setDT(bin_pops)

stmt <- "select * from doses_available"
doses_remaining <- dbGetQuery(db, stmt)
doses_remaining <- setDT(doses_remaining)
doses_remaining[, date_str := ymd(date_str)]

stmt <- "select * from doses_used"
doses_used <- dbGetQuery(db, stmt)
doses_used <- setDT(doses_used)
doses_used[, date_str := ymd(date_str)]
dbDisconnect(db)

doses = doses_remaining[, .(sim_day, date_str, dose, bin, std_remain=std_doses, urg_remain=urg_doses)][
  doses_used[, .(sim_day, date_str, dose, bin, std_used=std_doses, urg_used=urg_doses)], on = .(sim_day, date_str, dose, bin)][
    , .(sim_day, date_str, dose, bin, std_doses=(std_remain+std_used), urg_doses=(urg_remain+urg_used))]

bin_vax_delivery <- vax_hist[order(vax_date), .(n_doses = .N), by = .(vax_date, dose, p_age_bin)]
bin_vax_delivery[bin_pops, on = .(p_age_bin=bin_min), n_doses_p10k := n_doses * (1e4/bin_pop)]

tot_pop <- bin_pops[bin_min != 0, sum(bin_pop)]
tot_vax_delivery <- vax_hist[order(vax_date), .(n_doses = .N, cov = .N/tot_pop), by = .(vax_date, dose)]

dose_compare <- vax_hist[order(vax_date), .(doses_used=.N), by=.(date=vax_date, bin=p_age_bin, dose)][
  doses[order(date_str), .(doses_avail=sum(std_doses+urg_doses)), by=.(date=date_str, bin, dose)], on=.(date, bin, dose)]
dose_compare <- setnafill(x = dose_compare, fill = 0)
dose_compare[, diff := doses_avail - doses_used]
dose_compare[bin_pops, on = .(bin=bin_min), `:=` (
  doses_avail_p10k = doses_avail * (1e4/bin_pop),
  doses_used_p10k = doses_used * (1e4/bin_pop)
)]

# model vs. input total coverage by dose
tot_cov_by_dose_plot <- ggplot() +
  geom_line(data = exp.dt[order(exp, date, dose), .(tot_pop_vaxd = sum(tot_prop)), by = .(exp, date, dose)][
    ,cumsum := cumsum(tot_pop_vaxd), by=.(exp, dose)][
      ,.(date, location = sub("FL_like_(.*)", "\\1", exp), cumsum), by = .(exp, dose)][location == model_ref_loc],
    aes(x = date, y = cumsum), color = 'black') +
  geom_line(data = tot_vax_delivery[, .(vax_date, cumsum=cumsum(cov)), by=.(dose)][, dose := dose+1],
            aes(x = vax_date, y = cumsum), color = 'red') + 
  facet_grid(rows = vars(dose)) +
  ggtitle("Total coverage per dose (model vs. expectation)") +
  shr

# model vs. input total coverage by dose and bin
tot_cov_by_dose_bin_plot <- ggplot() +
  geom_line(data = prepended_doses.dt[, .(date, cumsum=cumsum(n_doses_p10k)), by=.(location, bin_min, dose)][location == model_ref_loc], aes(x = date, y = cumsum)) +
  geom_line(data = bin_vax_delivery[, .(vax_date, cumsum=cumsum(n_doses_p10k)), by=.(bin_min=p_age_bin, dose)][, dose := dose+1], aes(x = vax_date, y = cumsum), color = 'red') +
  facet_grid(rows=vars(dose), cols=vars(bin_min), labeller = labeller(bin_min = bin.labs)) +
  shr +
  xlab("Date (2020-2022)") + 
  ylab("Total doses delivered per 10K") +
  ggtitle("Cumulative administered doses (model vs. expectation)")

# model: cumul doses avail and used
cumul_doses_avail_v_used_plot <- ggplot(data = dose_compare[, .(doses_used=sum(doses_used), doses_avail=sum(doses_avail)), by=.(date)]) +
  geom_line(aes(x = date, y = cumsum(doses_avail))) +
  geom_line(aes(x = date, y = cumsum(doses_used)), linetype='dashed', color='red')

# model: diff between doses avail and used
doses_avail_v_used_diff_plot <- ggplot(data = dose_compare[, .(diff=sum(diff)), by=.(date)]) +
  geom_line(aes(x = date, y = diff))

# model: diff between doses avail and used by bin and dose 
doses_avail_v_used_diff_by_bin_dose_plot <- ggplot(data = dose_compare[bin!=0, .(date, bin, dose, diff)]) +
  geom_line(aes(x = date, y = diff)) +
  facet_grid(rows = vars(dose), cols = vars(bin), labeller = labeller(bin = bin.labs))

# doses avail from input file and model
doses_input_v_model_plot <- ggplot() +
  geom_line(data = prepended_doses.dt[, .(date, cumsum=cumsum(n_doses_p10k)), by=.(location, bin_min, dose)][location == model_ref_loc], aes(x = date, y = cumsum)) +
  geom_line(data = dose_compare[bin!=0, .(date, cumsum=cumsum(doses_avail_p10k)), by=.(bin_min=bin, dose)][, dose := dose+1], aes(x = date, y = cumsum), color='red', linetype='dashed') +
  facet_grid(rows=vars(dose), cols=vars(bin_min), labeller = labeller(bin_min = bin.labs))

ggsave(filename = paste0(figPath, "/tot_cov_by_dose.png"), plot = tot_cov_by_dose_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)
ggsave(filename = paste0(figPath, "/tot_cov_by_dose_bin.png"), plot = tot_cov_by_dose_bin_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)
ggsave(filename = paste0(figPath, "/cumul_doses_avail_v_used.png"), plot = cumul_doses_avail_v_used_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)
ggsave(filename = paste0(figPath, "/doses_avail_v_used_diff.png"), plot = doses_avail_v_used_diff_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)
ggsave(filename = paste0(figPath, "/doses_avail_v_used_diff_by_bin_dose.png"), plot = doses_avail_v_used_diff_by_bin_dose_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)
ggsave(filename = paste0(figPath, "/doses_input_v_model.png"), plot = doses_input_v_model_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)