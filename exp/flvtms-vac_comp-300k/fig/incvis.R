
.pkgs <- c("data.table", "ggplot2")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("digest.rds", "digest-key.rds", "vis_support.rda")),
  file.path("fig", "averted.png")
) else commandArgs(trailingOnly = TRUE)

#' comes key'd
eff.dt <- readRDS(.args[1])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, value
)]

scn.dt <- readRDS(.args[2])

load(.args[3])

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

filt <- if (interactive()) {
  expression(realization < 50)
} else expression(realization < 100)

p <- ggplot() + aes(
  x = date, y = value,
  color = active
) +
  geom_month_background(
    plt.dt[eval(filt)],
    font.size = 3,
    by = c(row = "outcome", col = "stockpile")
  ) +
  geom_spaghetti(
    mapping = aes(linetype = action, group = interaction(scenario, realization)),
    data = plt.dt[eval(filt)][!is.na(stockpile)]
  ) +
  geom_spaghetti(
    mapping = aes(linetype = action, group = interaction(scenario, realization)),
    data = plt.dt[eval(filt)][is.na(stockpile), .SD, .SDcol = -c("stockpile")]
  ) +
  facet_typical(scales = "free_y") +
  scale_y_incidence(
    name = "Per 10k, Incidence of ...", trans = "identity"
  ) +
  scale_x_null() +
  scale_color_scenario() +
  scale_linetype_scenario() +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 5*1.5, width = 2.25*4, bg = "white")
