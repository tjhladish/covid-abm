
.pkgs <- c(
  "data.table", "ggplot2", "scales",
  "ggh4x", "cabputils", "patchwork"
)

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("preference_wq.rds", "preference_nq.rds", "vocwindows.rds")),
  file.path("fig", "preference_ns.png")
))

load(.args[1])
wq_plot_dt <- readRDS(.args[2])[eval(seasfilter)]
nq_plot_dt <- readRDS(.args[3])[eval(seasfilter)]

scale_color_strategy <- rejig(
  scale_color_manual,
  name = "Vaccine Program",
  breaks = c("passive", "risk", "age", "ring"),
  labels = c(
    passive="Standard Program",
    ring="Ring Vaccination",
    age = "Age-Based Strategy",
    risk="Risk-Based Strategy"
  ),
  aesthetics = c("color", "fill"),
  values = c(passive = "black", ring = "#fb6502", risk = "#3b90db", age = "#209033"),
  guide = guide_legend(title.position = "top", title.hjust = 0.5)
)

mm.ref <- wq_plot_dt[,.(ymin = -25, ymax = 0), by = .(outcome)]
takeover.wins <- readRDS(.args[4])
tw <- takeover.wins[q == 0.5][CJ(measure, outcome = mm.ref$outcome), on=.(measure)][mm.ref, on=.(outcome)]
tw[, end := pmin(end, wq_plot_dt[, max(date)])]
tw[, mids := start + (end-start)/2 ]

plotter <- function(dt) ggplot(
  dt[nbadel != 0]
) + aes(
  date, neword, alpha = refv, fill = scenario
) + facet_grid(
  outcome ~ alloc, labeller = labeller(
    outcome = c(inf = "Infection", deaths = "Deaths")
  ), switch = "y"
) +
  geom_raster() +
  theme_minimal() + theme(
    legend.position = "bottom",
    legend.title = element_text(size = rel(0.65)),
    legend.text = element_text(size = rel(0.65)),
    legend.key.height = unit(0.1, "line"),
    panel.border = element_rect(fill = NA, linewidth = 0.25, color = "grey"),
    plot.margin = margin()
  ) +
  scale_alpha_continuous(guide = "none", range = c(0.2, 1)) +
  scale_y_continuous(name = NULL, breaks = NULL, expand = c(0, 0)) +
  scale_x_null(limits = dt[, {
    dr <- trunc.Date(range(date), "month")
    dr[2] <- trunc.Date(dr[2] + 31, "month") - 1
    dr
  }]) +
  scale_color_strategy(name = "Outcome Minimizing Program")

fs <- 3

voc.wins <- function (
  dt, al = 0.2, varcols = list(vocprev1 = "royalblue3",
  vocprev2 = "turquoise4", vocprev3 = "darkorchid3"),
  qs = c(0.5, 0.75, 0.95), ymin = 0, ymax = 1, laby = 0.05,
  vocs = c("α", "δ", "ο"), font.size = 12) {
  c(do.call(c, args = mapply(function(vc, tarq) geom_rect(aes(ymin = ymin,
    ymax = ymax, xmin = start, xmax = end), data = dt[measure == vc & q == tarq],
    alpha = al, inherit.aes = FALSE, fill = varcols[[vc]]
  ), vc = rep(names(varcols), each = length(qs)),
  tarq = rep(qs, times = length(varcols)), SIMPLIFY = FALSE)),
  if (length(vocs)) mapply(function(vc, lab) geom_text(aes(y = laby, x = mean(start) + 7),
    data = dt[measure == vc & q == 0.5], inherit.aes = FALSE,
    label = lab, color = varcols[[vc]],  hjust = 0, vjust = 0.5, size = font.size), vc = rep(names(varcols)),
    lab = vocs, SIMPLIFY = FALSE)
  )
}

p.bot <- ggplot() +
  geom_month_background(
    wq_plot_dt[, .(alloc, realization, date, value = nbadel) ],
    ymax = 1, ymin = 0, m.y = -0.1, y.y = -0.15, by = "alloc",
    font.size = fs*.75
  ) + facet_grid(. ~ alloc) + theme_minimal() + theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin()
  ) +
  voc.wins(
    tw, ymin = -0.5, ymax = 0, laby = -0.25, font.size = fs
  ) +
  coord_cartesian(ylim=c(-0.5, 0), expand = FALSE) +
  scale_x_null() + scale_y_continuous(NULL, guide = "none")

# want ordering by
#  - death outcomes
#  - the winning-est scenario descending
#  - within each present scenario, its advantage descending

# finout <- wq_plot_dt[
#   date == max(date) & outcome == "deaths",
#   .(alloc, realization, scenario, nbadel)
# ]
# finctr <- finout[, .N, by = .(alloc, scenario)][order(alloc, N)]
# finctr[, offset := c(0, head(cumsum(N), -1)), by = .(alloc)]
# jdt <- finout[finctr, on = .(alloc, scenario)]
# sdt <- jdt[, .SD[order(nbadel), .(realization, altord = seq_len(.N) + offset)], by = .(alloc, offset)]

reord <- function(dt) {
  finout <- dt[
    outcome == "deaths", .(
      nbadel = sum(nbadel), .N
    ), by = .(alloc, realization, scenario)
  ][, .SD[which.max(nbadel)], by = .(alloc, realization)]

  finctr <- finout[, .N, by = .(alloc, scenario)][order(alloc, N)]
  finctr[, offset := c(0, head(cumsum(N), -1)), by = .(alloc)]
  jdt <- finout[finctr[, .SD, .SDcols = -c("N")], on = .(alloc, scenario)]
  sdt <- jdt[, .SD[order(nbadel), .(realization, altord = seq_len(.N) + offset)], by = .(alloc, offset)]

  dt[sdt, on = .(alloc, realization), neword := altord]
  dt[]
}

wq_plot_dt <- reord(wq_plot_dt)
nq_plot_dt <- reord(nq_plot_dt)

p_wq <- (plotter(wq_plot_dt) / p.bot / guide_area()) + plot_layout(heights = c(1, .125, .05), guides = "collect")

# wq_plot_dt[date == max(date) & outcome == "deaths", altord := order(nbadel), by = .(alloc)]
# wq_plot_dt[, neword := altord[!is.na(altord)], by=.(alloc, realization)]

ggsave(
  gsub("\\.(.+)$", "_wq.\\1", tail(.args, 1)),
  p_wq, height = 5, width = 8, bg = "white"
)

p_nq <- (plotter(nq_plot_dt) / p.bot / guide_area()) + plot_layout(heights = c(1, .125, .05), guides = "collect")

ggsave(
  gsub("\\.(.+)$", "_nq.\\1", tail(.args, 1)),
  p_nq, height = 5, width = 8, bg = "white"
)
