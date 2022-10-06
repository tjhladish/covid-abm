#' SOURCES
#' Assume: COVID-19 vaccination in Sindh Province, Pakistan: a modelling study of health-impact and cost-effectiveness (SI)
#' https://doi.org/10.1371/journal.pmed.1003815
#' Empirical: UNICEF on 10 Aug 2022
#' https://www.unicef.org/supply/covid-19-vaccine-market-dashboard

.pkgs <- c("data.table")
stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumed to accessed via Rproj file, which puts at exp/active-vac
.args <- if (interactive()) c(
  "SES_stratified_doses.txt",
  "active_vax_counterfactual_doses.txt",
  "active_vac_doses"
) else commandArgs(trailingOnly = TRUE)

doses.in <- fread(.args[1])

dose_file_overwrite <- fread(.args[2])
dose_file_overwrite[, n_doses_p10k := 0]

cut_col <- function(col) {
  tmp <- copy(doses.in[dose_file_overwrite, on = .(date)][is_urg == 1 & bin_min == 5 & dose == 1, n_doses_p10k := col])
  return(tmp)
}

covax_hic = copy(doses.in[dose_file_overwrite, on = .(date)][is_urg == 1 & bin_min == 5 & dose == 1, n_doses_p10k := HIConly])
covax_mic = copy(doses.in[dose_file_overwrite, on = .(date)][is_urg == 1 & bin_min == 5 & dose == 1, n_doses_p10k := MIConly])
covax_lic = copy(doses.in[dose_file_overwrite, on = .(date)][is_urg == 1 & bin_min == 5 & dose == 1, n_doses_p10k := LIConly])

fwrite(covax_hic[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
       file = paste0(base::strsplit(tail(.args, 1), split = '\\.')[[1]][1], "_HIC_only.txt"), sep = " ")
fwrite(covax_mic[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
       file = paste0(base::strsplit(tail(.args, 1), split = '\\.')[[1]][1], "_MIC_only.txt"), sep = " ")
fwrite(covax_lic[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
       file = paste0(base::strsplit(tail(.args, 1), split = '\\.')[[1]][1], "_LIC_only.txt"), sep = " ")
