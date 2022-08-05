.pkgs <- c("data.table", "ggplot2", "patchwork")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes working directory at exp level,
#' though makefile will be invoked from w/in fig dir
.args = if (interactive()) c(
  "fig/validation.rds",
  "rcasedeath-florida.csv",
  "Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv",
  "COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv",
  "CDC_seroprev_long.csv",
  file.path("dose_data","trends_in_number_of_covid19_vaccinations_in_fl.csv"),
  file.path("fig", "vis_support.rda"),
  file.path("fig", "diagnosis.png")
) else commandArgs(trailingOnly=TRUE)

d <- readRDS(.args[1])

residue.dt <- melt(
  d[, .(realization, date, crcase, crdeath)],
  id.vars = c("realization", "date"), variable.name = "measure"
)[
  melt(ed[, .(date, crcase, crdeath)], id.vars = "date", variable.name = "measure"),
  on=.(date, measure)
][,
  residue := i.value - value
][,
  measure := gsub("^cr", "", measure)
][date < "2022-04-01"][!is.na(value)]

residue.dt$value <- residue.dt$i.value <- NULL
setnames(residue.dt, "residue", "value")

load(.args[7])

conserved <- list(
  scale_x_null(),
  scale_color_measure(),
  scale_shape_measure(),
  scale_alpha_measure(),
  theme_minimal()
)

p.core <- function(
  dt, ymax = NA, ymin = NA, ylog = FALSE, aesc = aes(color = measure),
  max.lines = 500, by=NULL
) ggplot(dt) +
  aes(date, value) + aesc +
  geom_month_background(dt, by = by, ymax = ymax, ymin = ymin, ylog = ylog) +
  geom_spaghetti(
    aes(y = value, group = interaction(realization, measure)),
    dt[!is.na(realization)], max.lines = max.lines
  )

p.cum.combo <- p.core(
  residue.dt,
  aesc = aes(color = measure), by = "measure"
) +
  facet_grid(measure ~ ., scales = "free_y") +
  scale_y_continuous(name = "Per 10k, Residual Cumulative Incidence of ...") +
  conserved +
  theme(legend.position = "none")

ggsave(tail(.args, 1), p.cum.combo, width = 12, height = 5, bg="white")
