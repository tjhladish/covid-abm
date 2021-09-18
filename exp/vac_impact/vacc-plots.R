
require(data.table)
require(ggplot2)

.args <- if (interactive()) c(
  "rawresults.rds",
  "vaccineimpact.png"
) else commandArgs(trailingOnly = TRUE)

mlt <- readRDS(.args[1])

# want the impact of introduction of variant, w/ and w/o vaccine

novac <- mlt[vac == 0][, -c("serial", "seed", "vac", "variable", "week")]
withvac <- mlt[vac == 1][, -c("serial", "seed", "vac", "variable", "week")]

var.impact <- withvac[
  date >= "2021-01-01" & var != "Rt"
][novac, on = .(realization, var, date, mutation), nomatch = 0]
var.impact[order(date), cvalue := cumsum(value), by=.(realization, var, mutation) ]
var.impact[order(date), i.cvalue := cumsum(i.value), by=.(realization, var, mutation) ]

# how much excess does variant cause?
var.impact[order(date), excess := value - i.value, by=.(realization, var, mutation) ]
var.impact[order(date), cexcess := cvalue - i.cvalue, by=.(realization, var, mutation) ]
# what is relative impact?
var.impact[order(date), ceff := -cexcess/i.cvalue, by=.(realization, var, mutation) ]

p <- ggplot(var.impact) + aes(
  date, ceff, color = mutation,
  group = interaction(realization, mutation)
) +
  facet_grid(var ~ ., scales = "free", labeller = labeller(
    var = c(c="Cases",d="Deaths")
  )) +
  geom_line(alpha = 0.1) +
  scale_x_date(
    name=NULL,
    date_breaks = "month", minor_breaks = NULL,
    date_labels = "%b %y"
  ) +
  scale_y_continuous("Relative Impact of Vaccine Introduction") +
  scale_color_continuous(
    "Variant", breaks = c(0,1), labels = c("None", "Some"),
    guide = guide_legend(override.aes = list(alpha=1))
  ) + theme_minimal() + theme(
    legend.position = c(1-0.05, 1-0.05),
    legend.justification = c(1, 1)
  ) + coord_cartesian(
    ylim = c(-1, NA)
  )

ggsave(tail(.args, 1), p, width = 7, height = 4, units = "in", dpi = 300, bg = "white")
