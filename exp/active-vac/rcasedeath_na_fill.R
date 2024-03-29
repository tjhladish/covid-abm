library(data.table)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/") }

.args <- if (interactive()) c(
  "rcasedeath-florida_copy.csv"
) else commandArgs(trailingOnly = TRUE)

d = fread(.args[1])
fwrite(x=d, file="./rcasedeath-florida_na_filled.csv", sep=",", na="NA", quote=F)
