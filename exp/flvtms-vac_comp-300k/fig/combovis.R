
.pkgs <- c("data.table", "ggplot2", "ggrepel", "patchwork")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("digest.rds", "digest-key.rds", "vis_support.rda")),
  file.path("fig", "combo-low.png")
) else commandArgs(trailingOnly = TRUE)

#' comes key'd
eff.dt <- readRDS(.args[1])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, value, averted, c.effectiveness
)]

scn.dt <- readRDS(.args[2])

load(.args[3])

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)[is.na(stockpile) | (stockpile == "low")]

rm(eff.dt)

aesbase <- aes(
  linetype = action,
  group = interaction(scenario, realization)
)

p.gen <- function(
  dt, yvar, aesbase, ytitle, scale_y,
  count = if (interactive()) 10 else 250,
  ymax = NA, show.end = FALSE
) {
  aesbase$y <-str2lang(yvar)
  loc.dt <- dt[realization < count][, c(key(dt), yvar), with = FALSE]
  return(ggplot() + aes_string(
    x = "date", y = yvar,
    color = "active"
  ) +
    geom_month_background(
      loc.dt,
      font.size = 3,
      by = c(row = "outcome"), value.col = yvar, ymax = ymax
    ) +
    geom_spaghetti(
      mapping = aesbase,
      data = loc.dt, #[!is.na(stockpile)],
      show.end = show.end
    ) +
    # geom_spaghetti(
    #   mapping = aesbase,
    #   data = loc.dt[is.na(stockpile), .SD, .SDcol = -c("stockpile")],
    #   show.end = show.end
    # ) +
    facet_typical(cols = NULL) +
    scale_y() +
    scale_x_null() +
    scale_color_scenario() +
    scale_linetype_scenario() +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = 8, margin = margin(b = 0.1, unit = "line")), # facet
      axis.text = element_text(size = 8), # ticks
      axis.title = element_text(size = 8, margin = margin(b = 0.1, unit = "line")) # axis label
    ))
}

p.eff <- p.gen(
  plt.dt[scenario != 1], "c.effectiveness", aesbase, "Test", scale_y_effectiveness,
  ymax = 1, show.end = TRUE
) + coord_cartesian(ylim = c(0, 1)) + theme(legend.position = "none")

p.ave <- p.gen(
  plt.dt[scenario != 1], "averted", aesbase, "Test", scale_y_averted
) + theme(legend.position = "none")

p.inc = p.gen(
  plt.dt, "value", aesbase, "Test", scale_y_incidence
) + theme(
  legend.position = "bottom",
  legend.margin = margin(),
  legend.spacing = unit(0.5, "line")
)

p.fin <- p.inc + p.ave + p.eff + guide_area() + plot_layout(
  design = "
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  DDD
", guides = 'collect')

ggsave(tail(.args, 1), p.fin, height = 6*1.5, width = 3.25*4, bg = "white")
