
.pkgs <- c("data.table", "scales", "ggplot2", "ggrepel", "patchwork", "cabputils", "geomtextpath")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

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
    "detection.rds",
    "vocpattern.rds",
    "vocwindows.rds"
  )),
  file.path("fig", "output", "everything_bw.png")
))

intfilter <- if (interactive()) expression(realization < 10) else expression(realization >= 0)

load(.args[1])
d <- readRDS(.args[2])[date <= vendday][eval(intfilter)]
ed <- readRDS(.args[3])[date <= vendday]
cdc = readRDS(.args[4])[date <= vendday]
hhsHosp = readRDS(.args[5])[date <= vendday][, .(date, hospInc) ]
seroprev = readRDS(.args[6])
vax = readRDS(.args[7])[date <= vendday]

ref.day0 <- d[, min(date)]

detect.dt <- readRDS(.args[8])[, date := day + ref.day0][date <= vendday]
voc.dt <- readRDS(.args[9])[date <= vendday]
takeover.win <- readRDS(.args[10])
takeover.win[, end := pmin(end, vendday)]
sd.dt <- d[realization == 1, .(date, value = sd, closed) ]

conserved <- list(
  scale_x_null(),
  scale_shape_measure(),
  theme_minimal(),
  theme(text = element_text(face = "bold")),
  scale_alpha_continuous(guide = "none", range = c(0.05, 1))
)

geom_grid <- function(
  ys, FUN = as.character
) {
  geom_texthline(
    aes(yintercept = ys, label = label),
    data = data.table(ys = ys, label = FUN(ys)), inherit.aes = FALSE,
    gap = TRUE, hjust = 0.325, color = "grey65", fontface = "bold"
  )
}

geom_month_end <- function(
    data,
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    m.labels = m.abb,
    font.size = 8, font.face = "bold",
    ytrans = "identity",
    datafn = tsref.dt,
    m.y = 0.0, y.y = 0.3,
    text.col = rep("darkgrey", length(col.cycle)),
    ...
) {
  names(text.col) <- c(tail(names(col.cycle), -1), names(col.cycle)[1])
  text.col[is.na(text.col)] <- "white"
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]
  dt[, col := text.col[mtype] ]

  xformer <- scales::as.trans(ytrans)

  xform <- function(ymax, ymin, dropto) {
    xformer$inverse(
      (xformer$transform(ymax)-xformer$transform(ymin))*dropto + xformer$transform(ymin)
    )
  }

  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (xformer$name %like% "log") 0 else -Inf, ymax = if (xformer$name %like% "prob") 1 else Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    ),
    geom_text(
      mapping = aes(x = mid, y = xform(ymax, ymin, m.y), label = m.abb[mon]),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, color = dt$col,
      size = font.size,
      vjust = if (y.y > 0) "bottom" else "top",
      fontface = font.face
    ),
    geom_text(
      mapping = aes(x = start+5, y = xform(ymax, ymin, y.y), label = yr),
      data = dt[yshow == TRUE],
      inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
      size = font.size, hjust = 0, fontface = font.face,
      vjust = if (y.y > 0) "top" else "bottom"
    )
  )
}

p.top <- ggplot() +
  geom_month_end(
    d[, .(realization, date, value = cov1) ],
    ymax = 1, ymin = 0, m.y = 0.05, y.y = 0.6
  ) + conserved + theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank()
  ) + coord_cartesian(ylim=c(0,0.6))

p.bot <- ggplot() +
  geom_month_end(
    d[, .(realization, date, value = cov1) ],
    ymax = 1, ymin = 0, m.y = -0.05, y.y = -0.6
  ) + conserved + theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank()
  ) + coord_cartesian(ylim=c(-0.6,0))

geom_month_bars <- function(
    data,
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    ytrans = "identity",
    datafn = tsref.dt,
    ...
) {
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]

  xformer <- scales::as.trans(ytrans)

  xform <- function(ymax, ymin, dropto) {
    xformer$inverse(
      (xformer$transform(ymax)-xformer$transform(ymin))*dropto + xformer$transform(ymin)
    )
  }

  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (xformer$name %like% "log") 0 else -Inf, ymax = if (xformer$name %like% "prob") 1 else Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    )
  )
}


p.core <- function(
  dt, ymax = NA, ymin = NA, ytrans = "identity", ins = list(), gridy
) ggplot(dt) +
  aes(date, value) +
  geom_grid(gridy) +
  geom_month_bars(
    dt, by = NULL, ymax = ymax, ymin = ymin, ytrans = ytrans
  ) +
  ins +
  stat_spaghetti(
    aes(sample = realization, alpha = after_stat(sampleN^-1)),
    data = \(d) d[!is.na(realization)],
    show.legend = TRUE
  )

sero.dt <- prepare(d[, .(realization, date, seroprev) ])
seroprev[, measure := "infection"]

geom_liner <- function(datafn) geom_text(
  aes(label = lbl, vjust = vj, hjust = hj, shape = NULL),
  data = datafn,
  fontface = "bold", size = 5
)

p.sero <- p.core(
  sero.dt[, measure := "infection"], ymin = 0.01, ymax = .99,
  gridy = c(0.25, 0.5, 0.75)
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
  scale_y_fraction() +
  conserved +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank()
  ) + coord_cartesian(clip = "off") +
  scale_alpha_continuous(guide = "none", range = c(0.05/7, 1))

inc.dt <- prepare(
  d[, .(realization, date, case = rcase, death = rdeath)],
  ed[, .(date, case = rcase, death = rdeath)]
)[date < "2022-04-01"][!is.na(value)]

p.inc.c <- p.core(
  inc.dt[measure == "case"], ymin = 1e-2, ymax = 1e2, ytrans = "log10",
  ins = voc.box(
    takeover.win, qs = c(0.5), vocs = c("\u03B1", "\u03B4", "\u03BF"),
    ymin = 0.01, ymax = 100, trans = "log10", laby = 0.85
  ),
  gridy = c(0.1, 1, 10)
) + geom_observation() +
  scale_y_incidence(trans = "log10", breaks = 10^((-2):2), labels = c("0.01", "0.1", "1", "10", "100")) +
  conserved +
  coord_cartesian(ylim = c(1e-2, 100), expand = FALSE, clip = "off") +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  scale_alpha_continuous(guide = "none", range = c(0.05/7, 1))

p.inc.d <- p.core(
  inc.dt[measure == "death"], ymin = 1e-2, ymax = 10, ytrans = "log10",
  ins = voc.box(
    takeover.win, qs = c(0.5), vocs = c("\u03B1", "\u03B4", "\u03BF"),
    ymin = 0.01, ymax = 10, trans = "log10", laby = 0.85
  ),
  gridy = c(0.1, 1)
) + geom_observation() +
  scale_y_incidence(trans = "log10", breaks = 10^((-2):1), labels = c("0.01", "0.1", "1", "10")) +
  conserved +
  coord_cartesian(ylim = c(1e-2, 10), expand = FALSE) +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

hinc.dt <- prepare(
  d[, .(realization, date = as.Date(date), hospInc) ], # vaxHosp, hospPrev, unvaxHosp,
  hhsHosp
)[date < "2022-04-01"][!is.na(value)]

p.hosp <- p.core(
  hinc.dt, ytrans = "log10", ymin = 1e-2, ymax = 3,
  ins = voc.box(
    takeover.win, qs = c(0.5), vocs = c("\u03B1", "\u03B4", "\u03BF"),
    ymin = 0.01, ymax = 3, trans = "log10", laby = 0.85
  ),
  gridy = c(0.03, 0.1, 0.3, 1)
) + geom_observation() +
  scale_y_incidence(trans = "log10", breaks = 10^sort(c((-2):0, log10(3)-(2:0))), labels = c("0.01", "0.03", "0.1", "0.3", "1", "3")) +#, )
  conserved +
  coord_cartesian(ylim = c(1e-2, 3), expand = FALSE) +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  scale_alpha_continuous(guide = "none", range = c(0.05/7, 1))

brk.dt <- rbind(prepare(
  d[, .(realization, date = as.Date(date), brkthru) ],
  cdc[, .(date = as.Date(date), brkthru = brkthruRatio) ]
)[date < vendday][!is.na(value)],
  data.table(realization = NA, date = as.Date("2020-02-01"), measure = "brkthru", value = NA)
)
brk.dt[!is.na(realization) & date < "2020-12-01", value := NA]

p.brk <- p.core(
  brk.dt, ymin = 0, ymax = 1,
  gridy = c(0.25, 0.5, 0.75)
) + geom_observation() +
  scale_y_fraction(name = "Fraction of Cases") +
  conserved +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) + coord_cartesian(clip = "off")

sd.dt <- d[realization == 1, .(date, value = sd, closed) ]

p.core <- function(
  dt, ymin = NA, ymax = NA, ytrans = "identity", gridy
) ggplot(dt) +
  aes(date, value) +
  geom_grid(gridy) +
  geom_month_bars(dt, ymin = ymin, ymax = ymax, ytrans = ytrans) +
  theme_minimal() + theme(text = element_text(face = "bold")) +
  scale_x_null()

schools <- data.table(
  start = as.Date(c("2020-02-01", "2020-09-01", "2021-08-10")),
  end = as.Date(c("2020-03-16", "2021-06-16", "2022-03-31")),
  level = c(0, .5, .2)
)

#' TODO school terms?
#' weekends?
p.sd <- p.core(
  sd.dt[, measure := "socialdist"], ymin = 0, ymax = 1,
  gridy = c(0.25, 0.5, 0.75)
) +
  geom_rect(
    aes(
      ymin = 0, ymax = 1, xmin = start, xmax = end
    ),
    data = function (dt) d[closed == 1, {
      spn = range(date)
      .(start = spn[1], end = spn[2])
    } ],
    inherit.aes = FALSE, alpha = 0.35, fill = "grey50"
  ) +
  geom_text(
    aes(
      y = 0.175, x = start+(end-start)/2, label = "Lockdown"
    ),
    data = function (dt) d[closed == 1, {
      spn = range(date)
      .(start = spn[1], end = spn[2])
    } ],
    inherit.aes = FALSE, angle = 90, hjust = 0, size = 6
  ) +
  # geom_rect(
  #   aes(
  #     ymin = level, ymax = 1, xmin = start - 0.5, xmax = end + 0.5
  #   ),
  #   data = schools,
  #   inherit.aes = FALSE, alpha = 0.25, fill = "dodgerblue"
  # ) +
  # geom_text(
  #   aes(
  #     y = 0.975, x = start+(end-start)/2, label = "School\nActivity"
  #   ),
  #   data = schools[2],
  #   inherit.aes = FALSE, vjust = 1, hjust = 0.5, size = 6
  # ) +
  geom_line() +
  scale_y_fraction(
    name = "Risk Threshold",
    sec.axis = sec_axis(
      name = "Per 10k,\nDaily Cases Reported",
      trans = function(x) 4*(x-0.5),
#      breaks = seq(0, 1, by=.25),
      labels = c("0.01", "0.1", "1", "10", "100")
    )
  ) + theme(panel.grid.major.y = element_blank()) +
  coord_cartesian(clip = "off")

seas.dt <- prepare(d[realization == 1, .(realization, date, seasonality) ])

scale_y_seasonality <- rejig(
  scale_y_continuous,
  name = "Seasonal Transmission Multiplier",
  limits = c(0.8, 1.2), breaks = seq(0.8, 1.2, by=.1),
  expand = expansion(
    mult = c(0, 0), add = c(0, 0)
  )
)

p.seas <- p.core(
  seas.dt, ymin = 0.8, ymax = 1.2, gridy = c(0.85, 1, 1.15)
) +
  geom_line() +
  scale_y_seasonality() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  coord_cartesian(clip = "off")

#' TODO make geom_month_background work for logit
p.voc <- p.core(
  voc.dt[!is.na(value)][value != 0], ymin = 0, ymax = 1,
  gridy = c(0.25, 0.5, 0.75)
) +
  voc.wins(takeover.win, qs = c(0.5, 0.75, 0.95), vocs = c("\u03B1", "\u03B4", "\u03BF")) +
  stat_spaghetti(aes(sample = realization, alpha = after_stat(sampleN^-1))) +
  scale_y_fraction(
    name = "Variant Fraction"
  ) + scale_alpha_continuous(range = c(0.05, 1)) +
  coord_cartesian(clip = "off")

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
  ymin = 0, ymax = 1,
  gridy=c(0.25, 0.5, 0.75)
) +
  aes(linetype = event, shape = event) +
  geom_line(data = function(dt) dt[realization == 1][value > 0]) +
  geom_text(mapping = aes(label=lab, linetype = NULL, shape = NULL), data = data.table(
    lab = c("Dose 1", "Dose 2", "Dose 3"),
    date = as.Date(c("2021-05-15", "2021-05-15", "2021-11-15")),
    value = c(.6, .15, .25),
    measure = "coverage"
  ), size = 6) +
  scale_y_fraction() +
  scale_linetype_manual(
    name = "Doses",
    values = c(cov1="dotted", cov2="dashed", cov3="solid"),
    labels = vlbls, guide = "none"
  ) +
  scale_shape_manual(
    name = "Doses",
    values = c(cov1=1, cov2=10, cov3=19),
    labels = vlbls,
    guide = "none"
  ) +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank()
  ) +
  coord_cartesian(clip = "off")

det.dt <- prepare(detect.dt[, .(
  realization = 1, date = day + ref.day0, asymp, mild, severe, crit
)])

setnames(det.dt, "measure", "outcome")
det.dt[, measure := "detection"]

p.detect <- p.core(
  det.dt, ymin = 0, ymax = 1, gridy = c(0.25, 0.5, 0.75)
) + aes(linetype = outcome) +
  geom_line() +
  scale_y_fraction(name = "P(Detect) Individual\nwith Outcome ...") +
  geom_text(mapping = aes(label=lab, linetype = NULL, shape = NULL, hjust = align), data = data.table(
    lab = c("Asymptomatic", "Mild", "Severe", "Critical"),
    date = as.Date(c("2020-07-15", "2020-05-10", "2020-04-15", "2020-08-22")),
    value = c(.18, .22, .65, .9),
    align = c(0.5, 1, 0.5, 0.5),
    measure = "detection"
  ), size = 6) +
  scale_linetype_manual(
    name = NULL, breaks = rev(c("asymp", "mild", "severe", "crit")),
    labels = c(asymp = "Asymptomic", mild = "Mild", severe = "Severe", crit = "Critical"),
    values = c(asymp = "dotted", mild = "dotdash", severe = "longdash", crit = "solid")
  ) +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank()
  ) + coord_cartesian(clip = "off")

p.res <- p.top + p.seas + p.detect + p.vax + p.sd +
  p.inc.c + p.hosp + p.inc.d + p.sero + p.brk +
  #p.cum.combo +
  p.bot +
  plot_layout(ncol = 1, heights = c(0.4, rep(1, 9), 0.4)) +
  plot_annotation(tag_levels = list(c(
    "", seas = "(a) Seasonality",
    det = "(b) Detection %", doses = "(c) Vaccine Coverage",
    sd = "(d) Societal Risk Perception",
    incc = "(e) Reported Cases (Daily, Log Scale)",
  #  cinc = "Per 10k, Cumulative Incidence of ...",
    hinc = "(f) Hospital Admissions (Daily, Log Scale)",
  incd = "(g) Excess Deaths (Weekly, Log Scale)",
    sero = "(h) Seroprevalence (Spike IgG)", brk = "(i) Breakthrough Infections in Vaccinees", ""
  ))) &
  theme(
    plot.tag.position = c(0, 0.97),
    axis.title = element_blank(), axis.text = element_blank(),
    plot.tag = element_text(size = rel(1.65), hjust = 0, vjust = 1)
  )

store(p.res, .args, width = 14, height = 20, bg = "white")
