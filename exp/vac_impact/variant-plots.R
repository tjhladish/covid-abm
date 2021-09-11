
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

var.impact <- variant[novariant, on = .(realization, var, date, vac)]
var.impact

ggplot(impact) + aes(date, value,
  color = linetype=c("NONE","VAC")[vac+1],
  group = interaction(serial,realization,vac,mutation)
) +
  facet_grid(var ~ ., scales = "free") +
  geom_line(alpha = 0.1) +
  theme_minimal()
