library(data.table)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/") }

.args = if (interactive()) c(
  "active_vac_raw_output.tgz",
  "ring_vax_counterfactual_dose_files"
) else commandArgs(trailingOnly = TRUE)

if (length(.args) != 2) {
  stop("Rscript active_plot_log_extract.R [path to zipped archive] [path to dir where dose files will be created]")
}

print("Extracting plot_log files from zipped archive...")
tmp_dir = paste0("./tmp", paste0(sample(c(sample(LETTERS, 10, replace = T), sample(0:9, 10, replace = T))), collapse = ''))
dir.create(tmp_dir)
untar(.args[1], exdir = tmp_dir)

print("Extracting data from plot_log files...")
dir.create(.args[2])
pb = txtProgressBar(min = 0, max = length(list.files(tmp_dir)), initial = 0, style = 3)
for (i in 1:length(list.files(tmp_dir))) {
  in_plot_log = tolower(list.files(tmp_dir)[i])
  in_serial = as.integer(sub("plot_log(.*).csv", "\\1", in_plot_log))
  
  if (!file.exists(file.path(tmp_dir, in_plot_log))) { next }
  
  active_cntfact <- fread("./state_based_counterfactual_doses.txt")
  passive_cntfact <- fread("./state_based_counterfactual_doses.txt")
  active_cntfact <- active_cntfact[ref_location == "FL"]
  passive_cntfact <- passive_cntfact[ref_location == "FL"]
  
  d <- fread(file.path(tmp_dir, in_plot_log))
  
  # for the active counterfactual, risk-based vax will be allocated the total number of doses ring vax but condensed over 30 days
  month_urg_deploy <- rep(d[, sum(urg_doses)] / 30, 30)
  month_urg_deploy <- c(month_urg_deploy, rep(0, active_cntfact[is_urg == 1 & date >= "2021-05-01" & bin_min == 5, .N] - length(month_urg_deploy)))
  active_cntfact[is_urg == 1 & date >= "2021-05-01" & bin_min == 5, n_doses_p10k := month_urg_deploy]
  
  # for the passive counterfactual, assign the daily number of urg doses to the passive campaign
  urg_doses <- d[, .(date, urg_doses)]
  tmp <- urg_doses[passive_cntfact[is_urg == 0 & dose == 1 & bin_min == 5], on = .(date)]
  passive_cntfact[is_urg == 1 & dose == 1 & bin_min == 5, n_doses_p10k := tmp[, urg_doses]]
  
  fwrite(active_cntfact, file = file.path(.args[2], paste0(in_serial, "_ring_vax_deployment_active_counterfactual_doses.txt")), sep = ' ')
  fwrite(passive_cntfact, file = file.path(.args[2], paste0(in_serial, "_ring_vax_deployment_passive_counterfactual_doses.txt")), sep = ' ')
  
  setTxtProgressBar(pb, i)
}
close(pb)
unlink(tmp_dir, recursive = T)
