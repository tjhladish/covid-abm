
.pkgs <- c("data.table", "ggplot2", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("vis_support.rda", "digest.rds", "digest-key.rds", "digest-ref.rds")),
  file.path("fig", "averted.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, averted
)]

scn.dt <- readRDS(.args[3])

plt.dt <- setkeyv(
  eff.dt[scn.dt[scenario != 1], on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

filt <- if (interactive()) {
  expression(realization < 50)
} else expression(realization < 100)

p <- ggplot(plt.dt) + aes(
  x = date, y = averted,
  color = quar, sample = realization
) +
  facet_grid(outcome ~ pas_vac + pas_alloc + act_vac + act_alloc) +
  geom_month_background(
    plt.dt,
    font.size = 3, value.col = "averted"
  ) +
  stat_spaghetti() +
  scale_y_averted() +
  scale_x_null() +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 5*1.5, width = 2.25*4, bg = "white")
