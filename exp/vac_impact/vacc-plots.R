
require(data.table)
require(ggplot2)
require(patchwork)

.args <- if (interactive()) c(
  "rawresults.rds",
  "vaccineimpact.png"
) else commandArgs(trailingOnly = TRUE)

mlt <- readRDS(.args[1])

# want the impact of introduction of variant, w/ and w/o vaccine

novac <- mlt[vac == 0][, -c("serial", "seed", "vac", "variable", "week")]
withvac <- mlt[vac == 1][, -c("serial", "seed", "vac", "variable", "week")]

vac.impact <- withvac[
  date >= "2021-01-01" & var != "Rt"
][novac, on = .(realization, var, date, mutation), nomatch = 0]

vac.impact[order(date), c(
  "cvalue", "i.cvalue", "excess", "cexcess", "ceff"
) := {
  cv <- cumsum(value)
  icv <- cumsum(i.value)
  .(cv, icv, value - i.value, cv-icv, -(cv-icv)/icv)
}, by=.(realization, var, mutation)][, vac := 1 ]

novariant <- mlt[mutation == 0][, -c("serial", "seed", "mutation", "variable", "week")]
variant <- mlt[mutation == 1][, -c("serial", "seed", "mutation", "variable", "week")]

var.impact <- variant[
  date >= "2021-01-01" & var != "Rt"
][novariant, on = .(realization, var, date, vac), nomatch = 0]

var.impact[order(date), c(
  "cvalue", "i.cvalue", "excess", "cexcess", "ceff"
) := {
  cv <- cumsum(value)
  icv <- cumsum(i.value)
  .(cv, icv, value - i.value, cv-icv, (cv-icv)/cv)
}, by=.(realization, var, vac)][, mutation := 1 ]

varscale <- function(show = TRUE) scale_linetype_manual(
  "VoCs",
  values = c("none" = "longdash", "present" = "solid"),
  guide = if (show) guide_legend(title.position = "top", override.aes = list(alpha=1)) else "none"
)

vacscale <- function(show = TRUE) scale_color_manual(
  "Vaccination",
  values = c("none" = "black", "present" = "dodgerblue"),
  guide = if (show) guide_legend(title.position = "top", override.aes = list(alpha=1)) else "none"
)

datescale <- function() scale_x_date(
  name=NULL,
  date_breaks = "month", minor_breaks = NULL,
  date_labels = "%b"
)

effplot <- function(
  fromdate = "2021-01-01",
  outcomes = c("c")
) ggplot(vac.impact[date >= fromdate][var %in% outcomes][, vac := 1 ]) + aes(
  date, ceff,
  linetype = c("none", "present")[mutation + 1],
  color = c("none", "present")[vac + 1],
  group = interaction(realization, mutation)
) +
  facet_grid(var ~ ., scales = "free", labeller = labeller(
    var = c(c="Cumulative Effectiveness on Symptomatic Infections,\nVaccination", d="Cumulative Effectiveness on Deaths,\nVaccination")
  ), switch = "y") +
  geom_line(linetype="solid", alpha = 0.025) +
  geom_line(
    aes(group = NULL),
    data = function(dt) dt[,
      .(ceff=median(ceff)),
      by=.(date, mutation, var, vac)
    ]
  ) +
  datescale() +
  scale_y_continuous(name = NULL) +
  vacscale(show = FALSE) + varscale() +
  coord_cartesian(
    ylim = c(NA, 1),
    expand = FALSE
  )

veffplot <- function(
  fromdate = "2021-01-01",
  outcomes = c("c")
) ggplot(var.impact[date >= fromdate][var %in% outcomes]) + aes(
  date, ceff,
  linetype = c("none", "present")[mutation + 1],
  color = c("none", "present")[vac + 1],
  group = interaction(realization, vac)
) +
  facet_grid(var ~ ., scales = "free", labeller = labeller(
    var = c(c="Cumulative Effectiveness on Symptomatic Infections,\nNo VoCs", d="Cumulative Effectiveness on Deaths,\nNo VoCs")
  ), switch = "y") +
  geom_line(linetype="solid", alpha = 0.025) +
  geom_line(
    aes(group = NULL),
    data = function(dt) dt[,
                           .(ceff=median(ceff)),
                           by=.(date, mutation, var, vac)
    ]
  ) +
  datescale() +
  scale_y_continuous(name = NULL) +
  vacscale() + varscale(show = FALSE) +
  coord_cartesian(
    ylim = c(NA, 1),
    expand = FALSE
  )

# legend.position = c(0+0.05, 1-0.05),
# legend.justification = c(0, 1),
# legend.direction = "horizontal"

epiplot <- function(
  fromdate = "2021-01-01",
  outcomes = c("c"),
  refpop = 375475/1e4,
  minrate = 1
) ggplot(mlt[date >= fromdate][var %in% outcomes]) + aes(
  date, value/refpop,
  color = c("none", "present")[vac + 1],
  linetype = c("none", "present")[mutation + 1],
  group = interaction(realization, mutation, vac)
) +
  facet_grid(var ~ ., scales = "free_y", labeller = labeller(
    var = c(c="Symptomatic Infections,\nper 10k", d="Deaths, per 10k")
  ), switch = "y") +
  geom_line(linetype="solid", alpha = 0.025) +
  geom_line(
    aes(group = NULL),
    data = function(dt) dt[,
      .(value=median(value)),
      by=.(date, vac, mutation, var)]
  ) +
  datescale() +
  scale_y_log10(
    name = NULL,
    labels = scales::label_number_si()
  ) +
  vacscale() +
  varscale() +
  coord_cartesian(
    ylim = c(minrate, NA),
    expand = FALSE
  )

outplot <- (
  epiplot() / effplot() #/ veffplot()
) + plot_layout(guides = "collect") & theme_minimal(
  base_size = 12 # base FONT size in pts
) & theme(
  legend.position = "bottom", strip.placement = "outside"
)

ggsave(tail(.args, 1), outplot, width = 11, height = 11, units = "in", dpi = 600, bg = "white")
