
.pkgs <- c("data.table", "ggplot2", "ggrepel")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("digest.rds", "digest-key.rds", "vis_support.rda")),
  file.path("fig", "effectiveness.png")
) else commandArgs(trailingOnly = TRUE)

#' comes key'd
#' TODO: identify proper reference point?
eff.dt <- readRDS(.args[1])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, c.effectiveness, c.averted, i.c.value
)]

eff.dt[,
  c0.effectiveness := c.averted/(i.c.value-i.c.value[1]),
  by = .(scenario, realization, outcome)
]

eff.dt$c.averted <- eff.dt$i.c.value <- NULL

scn.dt <- readRDS(.args[2])

load(.args[3])

plt.dt <- setkeyv(
  eff.dt[scn.dt[scenario != 1], on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

p.gen <- function(
  yvar = "c.effectiveness",
  yname = "Cumulative Effectiveness Against ...",
  aesspag = {
    tmp <- aes(
      linetype = action,
      group = interaction(scenario, realization)
    )
    tmp$y <- str2lang(yvar)
    tmp
  }
) ggplot() + aes_string(
  x = "date", y = yvar,
  color = "active"
) +
  geom_month_background(
    plt.dt,
    font.size = 3,
    by = c(row = "outcome", col = "stockpile"), value.col = yvar, ymax = 1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_spaghetti(
    mapping = aesspag,
    data = plt.dt[!is.na(stockpile)],
    yvar = yvar
    , show.end = TRUE
  ) +
  geom_spaghetti(
    mapping = aesspag,
    data = plt.dt[is.na(stockpile), .SD, .SDcol = -c("stockpile")],
    yvar = yvar
    , show.end = TRUE
  ) +
  facet_typical() +
  scale_y_effectiveness(name = yname) +
  scale_x_null() +
  scale_color_scenario() +
  scale_linetype_scenario() +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

peff <- p.gen()
ggsave(tail(.args, 1), peff, height = 5*1.5, width = 2.25*4, bg = "white")

peff0 <- p.gen(yvar = "c0.effectiveness", yname = "0-Baselined Cumulative Effectiveness Against ...")
ggsave(gsub("\\.","-b0.", tail(.args, 1)), peff0, height = 5*1.5, width = 2.25*4, bg = "white")
