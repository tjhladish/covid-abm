#!/usr/bin/env Rscript

.pkgs <- c("data.table", "MMWRweek", "ggplot2", "ggrepel", "patchwork")

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
  file.path("fig", "validation.png")
) else commandArgs(trailingOnly=TRUE)

d <- readRDS(.args[1])

#' doesn't seem to apply to any entries?
# is.na(d) = sapply(d, is.infinite)

#escambia_fraction    = 0.0153 # fraction of FL pop that lives in Escambia
tarpop <- gsub("^.+-(\\w+)\\.csv$", "\\1", .args[2])
pop <- c(florida = 21538187, escambia = 312212, dade = 2794464)[tarpop]
per10k <- 1e4/pop

ed <- fread(.args[2])
ed[, rcase := rcase * per10k ]
#' TODO fix upstream to only report on wday == 7? if so, need to change 7 below
ed[, rdeath := excess*per10k*7 ]
ed[wday(Date) != 7, rdeath := NA ]
ed[, crcase := cumsum(rcase) ]
ed[!is.na(rdeath), crdeath := cumsum(rdeath) ]
setnames(ed, "Date", "date")

cdc = fread(.args[3])[
  `Age group` == "all_ages_adj" &
  `Vaccine product` == "all_types" & outcome == "case"
][, date := {
  spl <- lapply(tstrsplit(`MMWR week`,"(?<=.{4})", perl = TRUE), as.integer)
  MMWRweek2Date(spl[[1]], spl[[2]])
}][,
   brkthruRatio := `Vaccinated with outcome`/(`Vaccinated with outcome` + `Unvaccinated with outcome`)
][,
  vaxOutcomeP10k := `Vaccinated with outcome` * (1e4/`Fully vaccinated population`)
]

hhsHosp = fread(.args[4])[
  state == 'FL'
][,
  date := as.Date(date)
][order(date)][,
  hospInc := previous_day_admission_adult_covid_confirmed * per10k
]

seroprev = fread(.args[5])[
  !is.na(est)
][,
  unique(.SD), keyby=date
][, run := {
  res <- rle(est)
  rep(1:length(res$lengths), res$lengths)
}][,.(
  start = date[1], end = date[.N],
  lower = lower[1], est = est[1], upper = upper[1]
), by=run
]

vax = dcast(melt(
  fread(.args[6], stringsAsFactors=F, header=T, skip=2)[order(Date), .(
  date = Date,
  dose1 = `Daily Count People Receiving Dose 1`,
  dose2 = `Daily Count of People Fully Vaccinated`,
  dose3 = `Daily Count People Receiving a First Booster Dose`
)][, paste0("cov",c(1,2,3)) := lapply(.(dose1, dose2, dose3), cumsum) ],
id.vars = "date")[, value := fifelse(variable %like% "dose", value * per10k, value/pop)],
  date ~ variable
)

load(.args[7])

conserved <- list(
  scale_x_null(),
  scale_color_measure(),
  scale_shape_measure(),
  scale_alpha_measure(),
  theme_minimal()
)

p.core <- function(dt, ymax, ymin, ylog = FALSE, aesc = aes(color = measure)) ggplot(dt) +
  aes(date, value) + aesc +
  geom_month_background(dt, by = NULL, ymax = ymax, ymin = ymin, ylog = ylog) +
  geom_spaghetti(
    aes(y = value, group = interaction(realization, measure)),
    dt[!is.na(realization)]
  )

sero.dt <- melt(
  d[, .(realization, date, seroprev) ],
  id.vars = c("realization", "date"),
  variable.name = "measure"
)

p.sero <- p.core(
  sero.dt, ymin = 0, ymax = 1,
  aesc = aes(color = "infection")
) + geom_crosshair(
  mapping = aes(
    x = start + (end+1-start)/2, xmin = start, xmax = end+1,
    ymax = upper, y = est, ymin = lower,
    shape = "observed"
  ),
  data = seroprev
) +
  geom_text(
    mapping = aes(label = lbl),
    data = function(dt) dt[
      date == "2021-12-01",
      .(measure = measure[1],
        date = date[1], value = mean(value)+.05*(fifelse(.BY == "cinf", 1, -1)),
        lbl = fifelse(.BY == "cinf", "Ever Infected", "Seropositive")
      ), by=.(as.character(measure))
    ], vjust = 0.5, hjust = 1
  ) +
  scale_y_fraction() +
  conserved +
  theme(
    legend.position = c(0+0.05, 1-0.15), legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
    legend.title = element_blank(),
    legend.margin = margin(0, 1, 1, 1),
    legend.key.height = unit(1.5, "line"),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm")
  )

inc.dt <- setkey(melt(rbind(
  d[, .(realization, date, rcase, rdeath)],
  ed[, .(realization = NA, date, rcase, rdeath)]
), id.vars = c("realization", "date"), variable.name = "measure")[
  ,
  measure := gsub("^r", "", measure)
], measure, realization, date)[date < "2022-04-01"][!is.na(value)]

p.inc <- p.core(
  inc.dt, ymin = 1e-2, ymax = 1e2, ylog = TRUE
) + geom_observation() +
  geom_text(
    mapping = aes(label = lbl),
    data = function(dt) dt[date == "2021-04-17",
                           .(date = date[1], value = mean(value)*(10^(fifelse(.BY == "case", 1, -1)/1.5)),
                             lbl = fifelse(.BY == "case", "Reported\nCases, Daily", "Excess\nDeaths, Weekly")),
                           by=.(measure)
    ], vjust = 0.5
  ) +
  scale_y_incidence(trans = "log10", breaks = 10^((-2):2), labels = c("0.01", "0.1", "1", "10", "100")) +
  conserved +
  coord_cartesian(ylim = c(1e-2, 100), expand = FALSE) +
  theme(legend.position = "none")

cum.dt <- setkey(melt(rbind(
  d[, .(realization, date, crcase, crdeath)],
  ed[, .(realization = NA, date, crcase, crdeath)]
), id.vars = c("realization", "date"), variable.name = "measure")[
  ,
  measure := gsub("^cr", "", measure)
], measure, realization, date)[date < "2022-04-01"][!is.na(value)]

p.cum.combo <- p.core(
  cum.dt, ymin = 1, ymax = 1e4, ylog = TRUE,
  aesc = aes(color = measure)
) + geom_observation() +
  geom_text(
    mapping = aes(label = lbl),
    data = function(dt) dt[date == "2021-04-17",
                           .(date = date[1], value = mean(value)*(10^(fifelse(.BY == "case", 1, -1)/1.75)),
                             lbl = fifelse(.BY == "case", "Cumulative Reported\nCases, Daily", "Cumulative Excess\nDeaths, Weekly")),
                           by=.(measure)
    ], vjust = 0.5
  ) +
  scale_y_log10(name = "Per 10k, Cumulative Incidence of ...") +
  conserved +
  coord_cartesian(ylim = c(1, 1e4), expand = FALSE) +
  theme(legend.position = "none")

hinc.dt <-setkey(melt(rbind(
  d[, .(realization, date = as.Date(date), vaxHosp, hospInc, hospPrev) ], # unvaxHosp,
  hhsHosp[, .(date, hospInc) ],
  fill = TRUE
), id.vars = c("realization", "date"), variable.name = "measure"),
  measure, realization, date)[date < "2022-04-01"][!is.na(value)]

p.hosp <- p.core(
  hinc.dt, ymin = 1e-2, ymax = 10, ylog = TRUE,
  aesc = aes(color = measure)
) + geom_observation() +
  geom_text(
    mapping = aes(label = lbl, hjust = ifelse(measure == "hospInc",1,0)),
    data = function(dt) dt[date == "2021-03-17",
                           .(measure = measure[1], date = date[1], value = max(mean(value), 1e-2)*(10^(fifelse(.BY == "hospInc", -1, 1)/3)),
                             lbl = c(hospPrev = "Hospital Occupancy,\nDaily", hospInc = "Total Admissions,\nDaily", vaxHosp = "Vaccinated Admissions,\nDaily")[unlist(.BY)]
                           ),
                           by=.(as.character(measure))
    ], vjust = 0.5
  ) +
  scale_y_log10(name = "Per 10k, ...") + # , breaks = 10^((-2):2), labels = c("0.01", "0.1", "1", "10", "100")
  conserved +
  coord_cartesian(ylim = c(1e-2, 10), expand = FALSE) +
  theme(legend.position = "none")

p.res <- p.sero + p.inc + p.cum.combo + p.hosp + plot_layout(nrow = 4)

ggsave(tail(.args, 1), p.res, width = 14, height = 11, bg = "white")

# TODO this seems like a more useful diagnostic view, but too complicated for lay readers
# residue.dt <- melt(
#   d[, .(realization, date, crcase, crdeath)],
#   id.vars = c("realization", "date"), variable.name = "measure"
# )[
#   melt(ed[, .(date, crcase, crdeath)], id.vars = "date", variable.name = "measure"),
#   on=.(date, measure)
# ][,
#   residue := i.value - value
# ][,
#   measure := gsub("^cr", "", measure)
# ][date < "2022-04-01"][!is.na(value)]
#
# residue.dt$value <- residue.dt$i.value <- NULL
# setnames(residue.dt, "residue", "value")

# p.cum <- function(dt, meas, ymax) ggplot(dt[measure == meas]) + aes(
#   date + 0.5, value, color = measure
# ) +
#   geom_month_background(dt[measure == meas], by = NULL, ymax = ymax) +
#   geom_spaghetti(
#     aes(y=value, group = interaction(realization, measure)),
#     dt[measure == meas][!is.na(realization)]
#   ) +
#   geom_observation() +
#   geom_text(
#     mapping = aes(label = lbl),
#     data = function(dt) dt[date == "2021-04-17",
#                            .(date = date[1], value = mean(value)*1.5,
#                              lbl = fifelse(.BY == "case", "Cumulative Reported\nCases, Daily", "Cumulative Excess\nDeaths, Weekly")),
#                            by=.(measure)
#     ], vjust = 0.5
#   ) +
#   scale_y_continuous(name = "Per 10k, Cumulative Incidence of ...", labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
#   conserved +
#   coord_cartesian(ylim = c(0, ymax), expand = FALSE) +
#   theme(legend.position = "none")
#
# p.cum.case <- p.cum(cum.dt, "case", 3000)
# p.cum.death <- p.cum(cum.dt, "death", 50)
