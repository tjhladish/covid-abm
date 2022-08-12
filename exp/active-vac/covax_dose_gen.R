#' SOURCES
#' Assume: COVID-19 vaccination in Sindh Province, Pakistan: a modelling study of health-impact and cost-effectiveness (SI)
#' https://doi.org/10.1371/journal.pmed.1003815
#' Empirical: UNICEF on 10 Aug 2022
#' https://www.unicef.org/supply/covid-19-vaccine-market-dashboard

.pkgs <- c("data.table", "wpp2019", "countrycode")
stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumed to accessed via Rproj file, which puts at exp/active-vac
.args <- if (interactive()) c(
  "covax_raw.csv",
  "active_vax_counterfactual_doses.txt",
  "covax_doses.txt"
) else commandArgs(trailingOnly = TRUE)

dose_file <- fread(.args[1])[,
  COVAX := as.integer(gsub(",","",COVAX))
] |> dcast(`Country/territory` + `mmm Year` ~ `WHO Preferred Vaccine Name`, value.var = "COVAX")

dropcols <- dose_file[, .SD |> lapply(\(col) all(is.na(col))), .SDcols = -c("Country/territory", "mmm Year")] |>
  unlist() |> which() |> names()

dose_file <- dose_file[, .SD, .SDcols = -c(dropcols)]

# assumes all individuals need two doses of all vaccines given
# despite Janssen being by label a 1-dose vaccine
dose_file[, doses := rowSums(.SD, na.rm = TRUE), .SDcols = -c("Country/territory", "mmm Year")]

data("popF"); data("popM")
refpop <- as.data.table(popF)[
  age != "0-4", .(pop1k = sum(`2020`)), keyby=.(country_code, name)
][as.data.table(popM)[
  age != "0-4", .(pop1k = sum(`2020`)), keyby=.(country_code, name)
]][, .(country_code, pop10k = (pop1k + i.pop1k)/10)]

covax.dt <- dose_file[, .(
  name = `Country/territory`,
  iso3 = countrycode(`Country/territory`, "country.name", "iso3n"),
  date = as.Date(paste0("01 ",`mmm Year`), "%d %b %Y"),
  doses
)][!(name %in% c("Humanitarian Buffer", "Kosovo"))]

globalpop10k <- refpop[country_code %in% unique(covax.dt$iso3), sum(pop10k)]
res.dt <- covax.dt[,
  .(dp10k = sum(doses)/globalpop10k, measure = "provisioned"), keyby=date
][CJ(
  date = seq(min(date), as.Date("2022-03-31"), by="day")
), .(
  date, dp10k = nafill(dp10k, "locf")
), on=.(date)][,.(
  date, provisioned = dp10k/with(rle(dp10k), rep(lengths, lengths))
)][which.max(provisioned > 0):.N]

# vs 10ks of 5+ used in Sindh model
sindh10k <- 4444.124
sindhcpd <- 4000 # courses per day, so doses x 2
sindhapprox <- data.table(
  date = seq(as.Date("2021-04-01"), length.out = 12, by="month"),
  assumed = c(rep(1, 3), rep(4, 3), rep(6, 3), rep(8, 3))*sindhcpd/sindh10k
)[CJ(
  date = seq(min(date), as.Date("2022-03-31"), by="day")
), on=.(date), .(date, assumedcourses = nafill(assumed, "locf")*2, assumeddoses = nafill(assumed, "locf"))]

res.dt[
  sindhapprox, on=.(date), c("assumedcourses","assumeddoses") := .(assumedcourses, assumeddoses)
][,
  assumedcourses := nafill(assumedcourses, fill = 0)][,
  assumeddoses := nafill(assumeddoses, fill = 0)
]

fwrite(res.dt, file = tail(.args, 1), sep = " ")

dose_file_overwrite <- fread(.args[2])
dose_file_overwrite[, n_doses_p10k := 0]
covax_assumed = copy(res.dt[dose_file_overwrite, on = .(date)][is_urg == 1 & bin_min == 5 & dose == 1, n_doses_p10k := assumedcourses])
covax_provisioned = copy(res.dt[dose_file_overwrite, on = .(date)][is_urg == 1 & bin_min == 5 & dose == 1, n_doses_p10k := provisioned])
fwrite(covax_assumed[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
       file = paste0(base::strsplit(tail(.args, 1), split = '\\.')[[1]][1], "_assumed.txt"), sep = " ")
fwrite(covax_provisioned[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)],
       file = paste0(base::strsplit(tail(.args, 1), split = '\\.')[[1]][1], "_provisioned.txt"), sep = " ")
