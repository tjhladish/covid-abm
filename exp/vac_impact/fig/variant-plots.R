
require(data.table)
require(ggplot2)

.args <- if (interactive()) c(
  "rawresults.rds",
  "variantimpact.png"
) else commandArgs(trailingOnly = TRUE)

mlt <- readRDS(.args[1])

# want the impact of introduction of variant, w/ and w/o vaccine

novariant <- mlt[mutation == 0][, -c("serial", "seed", "mutation", "variable", "week")]
variant <- mlt[mutation == 1][, -c("serial", "seed", "mutation", "variable", "week")]

var.impact <- variant[
  date >= "2021-01-01" & var != "Rt"
][novariant, on = .(realization, var, date, vac), nomatch = 0]
var.impact[order(date), cvalue := cumsum(value), by=.(realization, vac, var) ]
var.impact[order(date), i.cvalue := cumsum(i.value), by=.(realization, vac, var) ]

# how much excess does variant cause?
var.impact[order(date), excess := value - i.value, by=.(realization, vac, var) ]
var.impact[order(date), cexcess := cvalue - i.cvalue, by=.(realization, vac, var) ]
# what is relative impact?
var.impact[order(date), ceff := cexcess/i.cvalue, by=.(realization, vac, var) ]

p <- ggplot(var.impact) + aes(
  date, ceff, color = vac,
  group = interaction(realization, vac)
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
  scale_y_continuous("Relative Impact of Variant Introduction") +
  scale_color_continuous(
    "Vaccine", breaks = c(0,1), labels = c("None", "Some"),
    guide = guide_legend(override.aes = list(alpha=1))
  ) + theme_minimal() + theme(
    legend.position = c(0+0.05, 1-0.05),
    legend.justification = c(0, 1)
  )

ggsave(tail(.args, 1), p, width = 7, height = 4, units = "in", dpi = 300, bg = "white")
