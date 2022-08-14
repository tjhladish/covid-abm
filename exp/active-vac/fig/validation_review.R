#!/usr/bin/env Rscript

.pkgs <- c("data.table", "ggplot2", "ggrepel", "patchwork", "scales")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes working directory at exp level,
#' though makefile will be invoked from w/in fig dir
.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "fig/input/validation.rds",
  "fig/input/outcomes.rds",
  "fig/input/vaccines.rds",
  "fig/input/hospitals.rds",
  "fig/input/seroprev.rds",
  "fig/input/FLvaccines.rds",
  "fig/input/detection.rds",
  file.path("fig", "validation_review.png")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])
d <- readRDS(.args[2])[date <= endday]
ed <- readRDS(.args[3])
cdc = readRDS(.args[4])
hhsHosp = readRDS(.args[5])
seroprev = readRDS(.args[6])
vax = readRDS(.args[7])

conserved <- list(
  scale_x_null(),
  scale_color_measure(),
  scale_shape_measure(),
  scale_alpha_measure(guide = guide_legend(
    breaks = c("observed", "sample", "central"),
    labels = measlbls[c("observed", "sample", "central")],
    values = c(observed = 0.6, sample = 0.05, quantile = 0.5, central = 1)[c("observed", "sample", "central")],
    override.aes = list(
      linetype = c("blank", "solid", "solid"),
      size = c(10, 1, 3)
    )
  )),
  theme_minimal(),
  theme(text = element_text(face = "bold"))
)

p.core <- function(dt, ymax = NA, ymin = NA, ylog = FALSE, aesc = aes(color = measure)) ggplot(dt) +
  aes(date, value) + aesc +
  geom_month_background(dt, by = NULL, ymax = ymax, ymin = ymin, ylog = ylog) +
  geom_spaghetti(
    aes(y = value, group = interaction(measure, realization)),
    dt[!is.na(realization)]
  )

sero.dt <- prepare(d[, .(realization, date, seroprev) ])

geom_liner <- function(datafn) geom_text(
  aes(label = lbl, vjust = vj, hjust = hj),
  data = datafn,
  fontface = "bold", size = 5
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
  geom_liner(function(dt) dt[
    date == "2021-11-01",
    .(date = date[1], value = mean(value)+.05*(fifelse(.BY == "cinf", 1, -1)),
      lbl = fifelse(.BY == "cinf", "Ever Infected", "Seropositive"),
      vj = 1, hj = 0
    ), by = measure
  ]) +
  scale_y_fraction() +
  conserved +
#  scale_shape_discrete(guide = "none") +
  theme(
    legend.position = c(0+0.05, 1-0.175), legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
    legend.title = element_blank(),
    legend.margin = margin(0, 1, 1, 1),
    legend.key.height = unit(1.5, "line"),
    legend.key.width = unit(2, "line"),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(.5, "line"),
    legend.text = element_text(size = 16)
  )

inc.dt <- prepare(
  d[, .(realization, date, case = rcase, death = rdeath)],
  ed[, .(date, case = rcase, death = rdeath)]
)[date < "2022-04-01"][!is.na(value)]

p.inc <- p.core(
  inc.dt, ymin = 1e-2, ymax = 1e2, ylog = TRUE
) + geom_observation() +
  geom_liner(
    data = function(dt) dt[date == "2021-04-17", .(
      date = date[1], value = mean(value)*(10^(fifelse(.BY == "case", 1, -1)/1.5)),
      lbl = fifelse(.BY == "case", "Reported\nCases, Daily", "Excess\nDeaths, Weekly"),
      vj = 0.5, hj = 0.5
    ), by=.(measure)]
  ) +
  scale_y_incidence(trans = "log10", breaks = 10^((-2):2), labels = c("0.01", "0.1", "1", "10", "100")) +
  conserved +
  coord_cartesian(ylim = c(1e-2, 100), expand = FALSE) +
  theme(legend.position = "none")

cum.dt <- prepare(
  d[, .(realization, date, case = crcase, death = crdeath)],
  ed[, .(date, case = crcase, death = crdeath)]
)[date < "2022-04-01"][!is.na(value)]

p.cum.combo <- p.core(
  cum.dt, ymin = 1, ymax = 1e4, ylog = TRUE,
  aesc = aes(color = measure)
) + geom_observation() +
  geom_liner(
    function(dt) dt[date == "2021-04-17", .(
      date = date[1], value = mean(value)*(10^(-1/1.75)),
      lbl = fifelse(.BY == "case", "Cumulative Reported\nCases, Daily", "Cumulative Excess\nDeaths, Weekly"),
      vj = 0.5, hj = 0.5
    ), by=.(measure)
    ]
  ) +
  scale_y_log10(
    name = "Per 10k, Cumulative Incidence of ...",
    labels = number_format(scale_cut = cut_short_scale())) +
  conserved +
  coord_cartesian(ylim = c(1, 1e4), expand = FALSE) +
  theme(legend.position = "none")

hinc.dt <- prepare(
  d[, .(realization, date = as.Date(date), hospInc) ], # vaxHosp, hospPrev, unvaxHosp,
  hhsHosp[, .(date, hospInc) ]
)[date < "2022-04-01"][!is.na(value)]

p.hosp <- p.core(
  hinc.dt, ylog = TRUE, ymin = 1e-2, ymax = 3,
  aesc = aes(color = measure)
) + geom_observation() +
  geom_liner(function(dt) dt[date == "2021-03-17", .(
    measure = measure[1], date = date[1], value = max(mean(value), 1e-2)*(10^(fifelse(.BY == "hospInc", -1, 1)/3)),
    lbl = c(hospPrev = "Hospital Occupancy,\nDaily", hospInc = "Hospital Admissions,\nDaily", vaxHosp = "Vaccinated Admissions,\nDaily")[unlist(.BY)],
    vj = 1, hj = 0.5
    ), by=.(as.character(measure))
  ]) +
  scale_y_incidence(trans = "log10", breaks = 10^sort(c((-2):0, log10(3)-(2:0))), labels = c("0.01", "0.03", "0.1", "0.3", "1", "3")) +#, )
  conserved +
  coord_cartesian(ylim = c(1e-2, 3), expand = FALSE) +
  theme(legend.position = "none")

# min.breakthrough <- d[order(date),.SD[which.max(tot_std_doses + tot_urg_doses > 0)][1, date], by=realization][, min(V1)]

brk.dt <- prepare(
  d[, .(realization, date = as.Date(date), brkthru) ],
  cdc[, .(date = as.Date(date), brkthru = brkthruRatio) ]
)[date < "2022-04-01"][!is.na(value)]

p.brk <- p.core(
  brk.dt, ymin = 0, ymax = 1,
  aesc = aes(color = measure)
) + geom_observation() +
  geom_liner(function(dt) dt[date == "2021-08-01", .(
    date = date[1], value = mean(value)+0.15,
    lbl = c(brkthru = "Observed Fraction of Cases\nin Vaccinees, Weekly")[unlist(.BY)],
    vj = 0, hj = 1
  ), by=measure
  ]) +
  scale_y_fraction(name = "Fraction of Cases") +
  conserved +
  theme(legend.position = "none")

# plot(d$date, d$brkthruRatio, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,1), bty='n')
# shading()
# lines(cdc$date, cdc$brkthruRatio)
# lines(d$date, d$brkthruRatio, col='tan')
# lines(d$date, d$brkthruRatioAvg, col='tan4', lwd = 2)
# annotate('Breakthrough ratio')






p.res <- p.sero + p.inc + p.cum.combo + p.hosp + p.brk + plot_layout(nrow = 5)

ggsave(tail(.args, 1), p.res, width = 14, height = 13.5, bg = "white")

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
