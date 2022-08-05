
.pkgs <- c("data.table", "MMWRweek", "ggplot2", "ggrepel", "patchwork")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  "fig/validation.rds",
  "rcasedeath-florida.csv",
  "Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv",
  "COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv",
  "CDC_seroprev_long.csv",
  file.path("dose_data","trends_in_number_of_covid19_vaccinations_in_fl.csv"),
  "reporting_dump.csv",
  file.path("fig", "vis_support.rda"),
  file.path("fig", "modelinputs.png")
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

detect.dt <- fread(.args[7])

load(.args[8])

sd.dt <- d[realization == 1, .(date, value = sd, closed) ]

p.core <- function(
  dt, ymin = NA, ymax = NA, ylog = FALSE,
  aesc = aes(color = measure)
) ggplot(dt) +
  aes(date, value) + aesc +
  geom_month_background(dt, ymin = ymin, ymax = ymax, ylog = ylog)

ed.xform <- ed[, .(date, value = {
  tmp <- log10(rcase)
  # going to rescale to -2 to 2
  tmp[tmp < -2] <- NA
  tmp[tmp > 2] <- NA
  # now want 2=>0.5, -2=>0.5, then + 0.5
  tmp <- tmp/4 + 0.5
})]

#' TODO school terms?
#' weekends?
p.sd <- p.core(
  sd.dt, ymin = 0, ymax = 1,
  aesc = aes(color = "socialdist")
) +
  geom_point(data = ed.xform[date < "2022-04-01"], color = "black", alpha = 0.5) +
  geom_line() +
  geom_rect(
    aes(
      ymin = 0, ymax = 1, xmin = start, xmax = end,
      fill = "socialdist"
    ),
    data = function (dt) d[closed == 1, {
      spn = range(date)
      .(start = spn[1], end = spn[2])
    } ],
    inherit.aes = FALSE, alpha = 0.3
  ) + scale_y_fraction(
    name = "Risk Threshold",
    sec.axis = sec_axis(
      name = "Per 10k, Daily Cases Reported",
      trans = function(x) 4*(x-0.5),
#      breaks = seq(0, 1, by=.25),
      labels = c("0.01", "0.1", "1", "10", "100")
    )
  ) + theme_minimal() + scale_x_null() +
  scale_color_inputs()

seas.dt <- prepare(d[realization == 1, .(realization, date, seasonality) ])

scale_y_seasonality <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Seasonal Transmission Multiplier",
  limits = c(0.8, 1.2), breaks = seq(0.8, 1.2, by=.1),
  expand = expansion(
    mult = c(0, 0), add = c(0, 0)
  )
)

p.seas <- p.core(seas.dt, ymin = 0.8, ymax = 1.2) +
  geom_line() + theme_minimal() + scale_x_null() +
  scale_y_seasonality() +
  scale_color_inputs()

voc.dt <- prepare(
  d[, .(realization, date, vocprev1, vocprev2, vocprev3) ]
)

#' TODO make geom_month_background work for logit
p.voc <- p.core(voc.dt, ymin = 0, ymax = 1) +
  geom_spaghetti(
    aes(y = value, group = interaction(realization, measure)),
    voc.dt
  ) +
  scale_y_fraction(
    name = "Variant Fraction"
  ) +
  scale_alpha_measure(guide = "none") +
  scale_x_null() +
  scale_color_inputs() + theme_minimal()

vax.dt <- prepare(
  d[realization == 1, .(realization, date, cov1, cov2, cov3)],
  vax[, .(realization = 0, date, cov1, cov2, cov3)]
)

vlbls <- c(cov1="First", cov2="Second", cov3="Booster")

p.vax <- p.core(
  vax.dt,
  ymin = 0, ymax = 1, aesc = aes(color = "coverage")
) +
  aes(linetype = measure, shape = measure) +
  geom_point(data = function (dt) dt[realization == 0][value > 0], alpha = 0.2) +
  geom_line(data = function(dt) dt[realization == 1][value > 0]) +
  scale_color_inputs() +
  scale_x_null() +
  scale_y_fraction() +
  theme_minimal() +
  scale_linetype_manual(
    name = "Simulated Doses",
    values = c(cov1="dotted", cov2="dashed", cov3="solid"),
    labels = vlbls
  ) +
  scale_shape_manual(
    name = "Reported Doses",
    values = c(cov1=1, cov2=10, cov3=19),
    labels = vlbls,
    guide = guide_legend(order = 1)
  ) +
  theme(
    legend.position = c(0+0.05, 1-0.175), legend.justification = c(0, 1),
    legend.direction = "horizontal",
    legend.margin = margin(),
    legend.key.height = unit(1.5, "line"),
    legend.key.width = unit(1.5, "line"),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(.5, "line"),
  )

ref.day0 <- d[, min(date)]

det.dt <- prepare(detect.dt[, .(
  realization = 1, date = day + ref.day0, asymp, mild, severe, crit
)])

p.detect <- p.core(
  det.dt, ymin = 0, ymax = 1, aesc = aes(color = "detection")
) + aes(linetype = measure) +
  geom_line() +
  scale_x_null() +
  scale_y_fraction(name = "P(Detect) Individual with Outcome ...") +
  theme_minimal() +
  scale_linetype_manual(
    name = NULL, breaks = rev(c("asymp", "mild", "severe", "crit")),
    labels = c(asymp = "Asymptomic", mild = "Mild", severe = "Severe", crit = "Critical"),
    values = c(asymp = "dotted", mild = "dotdash", severe = "longdash", crit = "solid")
  ) +
  scale_color_manual(name = NULL, values = c("black"), guide = "none") +
  theme(
    legend.position = c(0.6, 0.6), legend.justification = c(0, 1)
  )

p.res <- p.sd + p.seas + p.voc + p.vax + p.detect + plot_layout(nrow = 5)

ggsave(tail(.args, 1), p.res, width = 14, height = 13.5, bg = "white")
