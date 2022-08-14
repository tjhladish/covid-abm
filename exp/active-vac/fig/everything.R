
.pkgs <- c("data.table", "scales", "ggplot2", "ggrepel", "patchwork")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "fig/input/validation.rds",
  "fig/input/outcomes.rds",
  "fig/input/vaccines.rds",
  "fig/input/hospitals.rds",
  "fig/input/seroprev.rds",
  "fig/input/FLvaccines.rds",
  "fig/input/detection.rds",
  file.path("fig", "output", "everything.png")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])
d <- readRDS(.args[2])[date <= endday]
ed <- readRDS(.args[3])
cdc = readRDS(.args[4])
hhsHosp = readRDS(.args[5])
seroprev = readRDS(.args[6])
vax = readRDS(.args[7])

ref.day0 <- d[, min(date)]

detect.dt <- readRDS(.args[8])[, date := day + ref.day0][date <= endday]

sd.dt <- d[realization == 1, .(date, value = sd, closed) ]

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

p.vis <- p.sero + p.inc + p.cum.combo + p.hosp + p.brk

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
p.voc <- p.core(
  voc.dt[!is.na(value)], ymin = 0, ymax = 1
) +
  geom_spaghetti(
    aes(y = value, group = interaction(measure, realization))
  ) +
  scale_y_fraction(
    name = "Variant Fraction"
  ) +
  scale_alpha_measure(guide = "none") +
  scale_x_null() +
  scale_color_inputs() + theme_minimal()

vax.mlt <- dcast(vax, date ~ dose, value.var = "cov")
setnames(vax.mlt, 2:4, paste0("cov", 1:3))
vax.dt <- prepare(
  d[realization == 1, .(realization, date, cov1, cov2, cov3)],
  vax.mlt[, .(realization = 0, date, cov1, cov2, cov3)]
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

p.res <- p.vis + p.sd + p.seas + p.voc + p.vax + p.detect + plot_layout(nrow = 10)

ggsave(tail(.args, 1), p.res, width = 14, height = 20, bg = "white")
