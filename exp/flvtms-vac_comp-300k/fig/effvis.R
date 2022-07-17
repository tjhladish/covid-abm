
.pkgs <- c("data.table", "ggplot2", "ggrepel")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("digest.rds", "digest-key.rds", "vis_support.rda")),
  file.path("fig", "effectiveness.png")
) else commandArgs(trailingOnly = TRUE)

#' comes key'd
eff.dt <- readRDS(.args[1])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, c.effectiveness
)]

scn.dt <- readRDS(.args[2])

load(.args[3])

plt.dt <- setkeyv(
  eff.dt[scn.dt[scenario != 1], on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

filt <- if (interactive()) {
  expression(realization < 10)
} else expression(realization < 750)

aesspag <- aes(
  y = c.effectiveness,
  linetype = action,
  group = interaction(scenario, realization)
)

p <- ggplot() + aes(
  x = date, y = c.effectiveness,
  color = active
) +
  geom_month_background(
    plt.dt[eval(filt)],
    font.size = 3,
    by = c(row = "outcome", col = "stockpile"), value.col = "c.effectiveness", ymax = 1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_spaghetti(
    mapping = aesspag,
    data = plt.dt[eval(filt)][!is.na(stockpile)]
    , show.end = TRUE
  ) +
  geom_spaghetti(
    mapping = aesspag,
    data = plt.dt[eval(filt)][is.na(stockpile), .SD, .SDcol = -c("stockpile")]
    , show.end = TRUE
  ) +
  facet_typical() +
  scale_y_effectiveness() +
  scale_x_null() +
  scale_color_scenario() +
  scale_linetype_scenario() +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 5*1.5, width = 2.25*4, bg = "white")
