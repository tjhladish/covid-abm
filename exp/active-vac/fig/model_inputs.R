
.pkgs <- c("data.table", "ggplot2", "ggrepel", "patchwork")

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
  file.path("fig", "model_inputs.png")
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

det.dt <- prepare(detect.dt[, .(
  realization = 1, date, asymp, mild, severe, crit
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
