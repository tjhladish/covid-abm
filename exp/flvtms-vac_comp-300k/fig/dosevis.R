
.pkgs <- c("data.table", "ggplot2", "ggrepel")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("digest-doses.rds", "digest-key.rds", "vis_support.rda")),
  file.path("fig", "doses.png")
) else commandArgs(trailingOnly = TRUE)

#' comes key'd
#' TODO: identify proper reference point?
doses.dt <- readRDS(.args[1])[
  date > "2020-12-01"
]

scn.dt <- readRDS(.args[2])

load(.args[3])

plt.dt <- setkeyv(
  doses.dt[scn.dt[!(action %in% c("none", "qonly"))], on=.(scenario)],
  union(key(doses.dt), colnames(scn.dt))
)

p <- ggplot(plt.dt) +
  aes(x = date, y = c.value, color = active) +
  geom_month_background(
    plt.dt,
    font.size = 3,
    by = c("stockpile"),
    value.col = "c.value"
  ) +
  geom_spaghetti(
    mapping = aes(
      y = c.value,
      linetype = action,
      group = interaction(scenario, realization)
    ),
    data = plt.dt
  ) +
  facet_typical(rows = NULL) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_cdoses() +
  scale_x_null() +
  scale_color_scenario() +
  scale_linetype_scenario() +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )


ggsave(tail(.args, 1), p, height = 5*1.5, width = 2.25*4, bg = "white")
