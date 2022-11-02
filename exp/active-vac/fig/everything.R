
.pkgs <- c("data.table", "scales", "ggplot2", "ggrepel", "patchwork", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = commandArgs(
  trailingOnly = TRUE, c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c(
    "validation.rds",
    "outcomes.rds",
    "vaccines.rds",
    "hospitals.rds",
    "seroprev.rds",
    "FLvaccines.rds",
    "detection.rds"
  )),
  file.path("fig", "output", "everything.png")
))

load(.args[1])
d <- readRDS(.args[2])[date <= vendday]
ed <- readRDS(.args[3])[date <= vendday]
cdc = readRDS(.args[4])[date <= vendday]
hhsHosp = readRDS(.args[5])[date <= vendday][, .(date, hospInc) ]
seroprev = readRDS(.args[6])
vax = readRDS(.args[7])[date <= vendday]

ref.day0 <- d[, min(date)]

detect.dt <- readRDS(.args[8])[, date := day + ref.day0][date <= vendday]

sd.dt <- d[realization == 1, .(date, value = sd, closed) ]

conserved <- list(
  scale_x_null(),
  scale_color_measure(),
  scale_shape_measure(),
  theme_minimal(),
  theme(text = element_text(face = "bold"))
)

p.core <- function(dt, ymax = NA, ymin = NA, ytrans = "identity") ggplot(dt) +
  aes(date, value, color = measure) +
  geom_month_background(
    dt, by = NULL, ymax = ymax, ymin = ymin, ytrans = ytrans
  ) +
  stat_spaghetti(
    aes(sample = realization),
    data = \(d) d[!is.na(realization)],
    show.legend = TRUE
  )

sero.dt <- prepare(d[, .(realization, date, seroprev) ])

geom_liner <- function(datafn) geom_text(
  aes(label = lbl, vjust = vj, hjust = hj, shape = NULL),
  data = datafn,
  fontface = "bold", size = 5
)

seroprev[, measure := "infection"]

p.sero <- p.core(
  sero.dt[, measure := "infection"], ymin = 0.01, ymax = .99
) + aes(shape = after_stat(spaghetti)) + geom_crosshair(
  mapping = aes(
    x = start + (end+1-start)/2, xmin = start, xmax = end+1,
    ymax = upper, y = est, ymin = lower
  ),
  data = seroprev
) + geom_observation(aes(
    x = start + (end+1-start)/2,
    y = est, shape="observation"
  ), data = seroprev) +
  geom_liner(function(dt) dt[
    date == "2021-11-01",
    .(date = date[1], value = mean(value)+.05*(fifelse(.BY == "cinf", 1, -1)),
      lbl = fifelse(.BY == "cinf", "Ever Infected", "Seropositive"),
      vj = 1, hj = 0
    ), by = measure
  ]) +
  scale_y_fraction(
#    breaks = c(0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99), limits = c(0.01, 0.99)
  ) +
  conserved +
#  coord_trans(y="logit") +
  scale_alpha(guide = "none") +
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
  inc.dt, ymin = 1e-2, ymax = 1e2, ytrans = "log10"
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
  cum.dt, ymin = 1, ymax = 1e4, ytrans = "log10"
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
  hhsHosp
)[date < "2022-04-01"][!is.na(value)]

p.hosp <- p.core(
  hinc.dt, ytrans = "log10", ymin = 1e-2, ymax = 3
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
  theme(legend.position = "none") +
  scale_alpha_continuous(range = c(0.01, 1))

# min.breakthrough <- d[order(date),.SD[which.max(tot_std_doses + tot_urg_doses > 0)][1, date], by=realization][, min(V1)]

brk.dt <- prepare(
  d[, .(realization, date = as.Date(date), brkthru) ],
  cdc[, .(date = as.Date(date), brkthru = brkthruRatio) ]
)[date < "2022-04-01"][!is.na(value)]

p.brk <- p.core(
  brk.dt, ymin = 0, ymax = 1
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
  dt, ymin = NA, ymax = NA, ytrans = "identity"
) ggplot(dt) +
  aes(date, value, color = measure) +
  geom_month_background(dt, ymin = ymin, ymax = ymax, ytrans = ytrans) +
  theme_minimal() + theme(text = element_text(face = "bold")) +
  scale_x_null()

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
  sd.dt[, measure := "socialdist"], ymin = 0, ymax = 1
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
  ) + scale_color_inputs()

seas.dt <- prepare(d[realization == 1, .(realization, date, seasonality) ])

scale_y_seasonality <- rejig(
  scale_y_continuous,
  name = "Seasonal Transmission Multiplier",
  limits = c(0.8, 1.2), breaks = seq(0.8, 1.2, by=.1),
  expand = expansion(
    mult = c(0, 0), add = c(0, 0)
  )
)

p.seas <- p.core(seas.dt, ymin = 0.8, ymax = 1.2) +
  geom_line() +
  scale_y_seasonality() +
  scale_color_inputs()

voc.dt <- prepare(
  d[, .(realization, date, vocprev1, vocprev2, vocprev3) ]
)

#' TODO make geom_month_background work for logit
p.voc <- p.core(
  voc.dt[!is.na(value)], ymin = 0, ymax = 1
) + stat_spaghetti(aes(sample = realization)) +
  scale_y_fraction(
    name = "Variant Fraction"
  ) +
  scale_color_inputs()

vax.mlt <- dcast(vax, date ~ dose, value.var = "cov")
setnames(vax.mlt, 2:4, paste0("cov", 1:3))
vax.dt <- prepare(
  d[realization == 1, .(realization, date, cov1, cov2, cov3)],
  vax.mlt[, .(realization = 0, date, cov1, cov2, cov3)]
)

vlbls <- c(cov1="First", cov2="Second", cov3="Booster")

setnames(vax.dt, "measure", "event")
vax.dt[, measure := "coverage"]
p.vax <- p.core(
  vax.dt,
  ymin = 0, ymax = 1
) +
  aes(linetype = event, shape = event) +
  geom_point(data = function (dt) dt[realization == 0][value > 0], alpha = 0.2) +
  geom_line(data = function(dt) dt[realization == 1][value > 0]) +
  scale_color_inputs() +
  scale_y_fraction() +
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

det.dt <- prepare(detect.dt[, .(
  realization = 1, date = day + ref.day0, asymp, mild, severe, crit
)])

setnames(det.dt, "measure", "outcome")
det.dt[, measure := "detection"]

p.detect <- p.core(
  det.dt, ymin = 0, ymax = 1
) + aes(linetype = outcome) +
  geom_line() +
  scale_y_fraction(name = "P(Detect) Individual with Outcome ...") +
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

store(.args, p.res, width = 14, height = 20, bg = "white")
