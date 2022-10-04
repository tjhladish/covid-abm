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

cut_col <- function(urg, col) {
  tmp <- copy(doses.in[dose_file_overwrite, on = .(date)][is_urg == urg & bin_min == 5 & dose == 1, n_doses_p10k := get(col)])
  return(tmp)
}

cut_col(1, "HIConly")[is_urg == 1 & bin_min == 5 & dose == 1]

covax_hic_std = cut_col(0, "HIConly")
covax_hic_urg = cut_col(1, "HIConly")

covax_mic_std = cut_col(0, "MIConly")
covax_mic_urg = cut_col(1, "MIConly")

covax_lic_std = cut_col(0, "LIConly")
covax_lic_urg = cut_col(1, "LIConly")

write_out <- function(dt, suffix) {
  fwrite(dt[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
         file = paste0(base::strsplit(tail(.args, 1), split = '\\.')[[1]][1], suffix), sep = " ")
}

write_out(covax_hic_std, "_HIC_std.txt")
write_out(covax_hic_urg, "_HIC_urg.txt")

write_out(covax_mic_std, "_MIC_std.txt")
write_out(covax_mic_urg, "_MIC_urg.txt")

write_out(covax_lic_std, "_LIC_std.txt")
write_out(covax_lic_urg, "_LIC_urg.txt")