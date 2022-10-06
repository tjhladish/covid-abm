#' SOURCES
#' Assume: COVID-19 vaccination in Sindh Province, Pakistan: a modelling study of health-impact and cost-effectiveness (SI)
#' https://doi.org/10.1371/journal.pmed.1003815
#' Empirical: UNICEF on 10 Aug 2022
#' https://www.unicef.org/supply/covid-19-vaccine-market-dashboard

.pkgs <- c("data.table")
stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumed to accessed via Rproj file, which puts at exp/active-vac
.args <- if (interactive()) c(
  "SES_stratified_doses_v2.txt",
  "active_vax_counterfactual_doses.txt"
) else commandArgs(trailingOnly = TRUE)

doses.in <- fread(.args[1])

dose_file_overwrite <- fread(.args[2])
dose_file_overwrite[, n_doses_p10k := 0]

cut_col <- function(urg, col) {
  tmp <- copy(doses.in[dose_file_overwrite, on = .(date)][is_urg == urg & bin_min == 5 & dose == 1, n_doses_p10k := get(col)])
  return(tmp)
}

covax_usa_std = cut_col(0, "US")
covax_usa_urg = cut_col(1, "US")

covax_hic_std = cut_col(0, "HIC")
covax_hic_urg = cut_col(1, "HIC")

covax_mic_std = cut_col(0, "MIC")
covax_mic_urg = cut_col(1, "MIC")

covax_lic_std = cut_col(0, "LIC")
covax_lic_urg = cut_col(1, "LIC")

write_out <- function(dt, suffix) {
  fwrite(dt[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
         file = paste0("active_vac_doses", suffix), sep = " ")
}

write_out(covax_usa_std, "_USA_std.txt")
write_out(covax_usa_urg, "_USA_urg.txt")

write_out(covax_hic_std, "_HIC_std.txt")
write_out(covax_hic_urg, "_HIC_urg.txt")

write_out(covax_mic_std, "_MIC_std.txt")
write_out(covax_mic_urg, "_MIC_urg.txt")

write_out(covax_lic_std, "_LIC_std.txt")
write_out(covax_lic_urg, "_LIC_urg.txt")
