#' SOURCES
#' Assume: COVID-19 vaccination in Sindh Province, Pakistan: a modelling study of health-impact and cost-effectiveness (SI)
#' https://doi.org/10.1371/journal.pmed.1003815
#' Empirical: UNICEF on 10 Aug 2022
#' https://www.unicef.org/supply/covid-19-vaccine-market-dashboard

.pkgs <- c("data.table", "wpp2019", "countrycode")
stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumed to accessed via Rproj file, which puts at exp/active-vac
.args <- if (interactive()) c(
  "unicef_raw.csv",
  "wb_inc.csv",
  "covax_doses.txt"
) else commandArgs(trailingOnly = TRUE)

dose_file <- fread(.args[1])[,
  c("COVAX", "bilat", "other") := .(
    as.integer(gsub(",","",COVAX)),
    as.integer(gsub(",","",`Bilateral/multilateral agreements`)),
    as.integer(gsub(",","",`Donations`)) + as.integer(gsub(",","",`Unknown`)) + as.integer(gsub(",","",`AVAT`))
  )
] |> melt.data.table(measure.vars = c("COVAX", "bilat", "other")) |> dcast(
  `Country/territory` + `mmm Year` + variable ~ `WHO Preferred Vaccine Name`,
  value.var = "value"
)

dropcols <- dose_file[, .SD |> lapply(\(col) all(is.na(col))), .SDcols = -c("Country/territory", "mmm Year")] |>
  unlist() |> which() |> names()

dose_file <- dose_file[, .SD, .SDcols = -c(dropcols)]

# assumes all individuals need two doses of all vaccines given
# despite e.g. Janssen being by label a 1-dose vaccine
dose_file <- dose_file[,
  doses := rowSums(.SD, na.rm = TRUE), .SDcols = -c("Country/territory", "mmm Year", "variable")
][!(`Country/territory` %in% c("Kosovo", "Humanitarian Buffer", "Hong Kong SAR and Macao")), .(
  name = `Country/territory`,
  date = as.Date.character(`mmm Year`, format = "%m/%d/%Y"),
  variable,
  iso3 = countrycode(`Country/territory`, "country.name", "iso3n"),
  doses
)]

data("popF"); data("popM")
refpop <- as.data.table(popF)[
  age != "0-4", .(pop1k = sum(`2020`)), keyby=.(country_code, name)
][as.data.table(popM)[
  age != "0-4", .(pop1k = sum(`2020`)), keyby=.(country_code, name)
]][, .(country_code, pop10k = (pop1k + i.pop1k)/10)]

wb.dt <- fread(.args[2])[!(country %in% c("Kosovo", "Channel Islands")), iso3 := countrycode(`country`, "country.name", "iso3n")]
wb.dt[country == "Türkiye", iso3 := 792]

scens <- list(
  US = expression(country == "United States"),
  HIC = expression(category == "HIC"),
  LIC = expression(category == "LIC"),
  UMIC = expression(category == "UMIC"),
  LMIC = expression(category == "LMIC"),
  MIC = expression(!(category %in% c("HIC", "LIC")))
)

covax.only <- dose_file[variable != "COVAX", sum(doses), by=.(iso3)][V1 == 0, iso3]

isos <- c(
  lapply(scens, function(e) wb.dt[eval(e), iso3]),
  list(COVAX = covax.only)
)

pop10k <- lapply(
  isos, function(is) refpop[country_code %in% is, sum(pop10k)]
)

globalpop10k <- refpop[country_code %in% unique(dose_file$iso3), sum(pop10k)]

doses.dt <- rbindlist(
  mapply(function(isogrp, pop) dose_file[
    iso3 %in% isogrp, .(dp10k = sum(doses)/pop), keyby=.(date)
  ], isogrp = isos, pop = pop10k, SIMPLIFY = FALSE),
  idcol = "grp"
)

covax.dt <- doses.dt[CJ(
  date = seq(min(date, na.rm = TRUE), as.Date("2022-03-31"), by="day"),
  grp = unique(grp)
), .(
  date, grp, dp10k
), on=.(date, grp)][order(date), dp10k := nafill(
  nafill(dp10k, "locf"), "nocb"
), by=.(grp) ][,.(
  date, provisioned = dp10k/with(rle(dp10k), rep(lengths, lengths))
), by=grp] |> dcast(date ~ grp, value.var = "provisioned")

# vs 10ks of 5+ used in Sindh model
sindh10k <- 4444.124
sindhcpd <- 4000 # courses per day, so doses x 2
sindhapprox <- data.table(
  date = seq(as.Date("2021-04-01"), length.out = 12, by="month"),
  assumed = c(rep(1, 3), rep(4, 3), rep(6, 3), rep(8, 3))*sindhcpd/sindh10k
)[CJ(
  date = seq(min(date), as.Date("2022-03-31"), by="day")
), on=.(date), .(date, assumedcourses = nafill(assumed, "locf")*2, assumeddoses = nafill(assumed, "locf"))]

covax.dt[
  sindhapprox, on=.(date), c("assumedcourses","assumeddoses") := .(assumedcourses, assumeddoses)
][,
  assumedcourses := nafill(assumedcourses, fill = 0)][,
  assumeddoses := nafill(assumeddoses, fill = 0)
]

fwrite(covax.dt, file = tail(.args, 1), sep = " ")

#' @examples
#' require(ggplot2)
#' ggplot(melt(covax.dt, id.var = "date")[, cvalue := cumsum(value), by=variable]) +
#'  aes(date, cvalue, color = variable) +
#'  geom_line(data = \(dt) dt[variable %in% c("US","HIC","MIC","LIC", "COVAX")]) +
#'  theme_minimal() + theme(legend.pos = c(0, 1), legend.jus = c(0, 1)) +
#'  scale_x_date()
