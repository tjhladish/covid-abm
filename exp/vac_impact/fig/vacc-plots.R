
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
  .(cv, icv, value - i.value, icv-cv, (icv-cv)/icv)
}, by=.(realization, var, mutation)][, 
  vac := 1
][,
  outcome := c(c="Symptomatic Infections", d="Deaths")[var]
]

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
}, by=.(realization, var, vac)][, 
  mutation := 1
][,
  outcome := c(c="Symptomatic Infections", d="Deaths")[var]
]

mlt[,
    outcome := c(c="Symptomatic Infections", d="Deaths")[var]
]

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
  breaks = function(l) {
    rl <- round(as.POSIXlt(l), "month")
    as.Date(seq(rl[1] + 60*60*24*15, rl[2], by="month"))
  },
  labels = function(bs) {
    gsub("^(.).+$","\\1",month.abb[month(bs)])
  },
  minor_breaks = NULL
)

datebg <- function(ylims = c(0,Inf), bandcolor = "grey") list(
  geom_rect(
    aes(xmin = start, xmax=end, ymax=ylims[2], ymin=ylims[1]),
    data = function(dt) {
      dr <- round(as.POSIXlt(dt[, range(date)]), "month")
      pts <- seq(dr[1], dr[2], by="month")
      if (length(pts) %% 2) pts <- head(pts, -1)
      data.table(
        start = as.Date(pts[seq(1, length(pts), by=2)]),
        end = as.Date(pts[seq(2, length(pts), by=2)])
      )
    },
    fill = bandcolor, alpha = 0.1,
    inherit.aes = FALSE, show.legend = FALSE
  ),
  datescale()
)

effplot <- function(
  fromdate = "2021-01-01",
  outcomes = c("c", "d")
) ggplot(vac.impact[date >= fromdate][var %in% outcomes]) + aes(
  date, ceff,
  linetype = c("none", "present")[mutation + 1],
  color = c("none", "present")[vac + 1],
  group = interaction(realization, mutation)
) +
  facet_grid(. ~ outcome) +
  datebg(ylims = c(-Inf,Inf)) +
  geom_line(linetype="solid", alpha = 0.025) +
  geom_line(
    aes(group = NULL),
    data = function(dt) dt[,
      .(ceff=median(ceff)),
      by=.(date, mutation, outcome, vac)
    ]
  ) +
  scale_y_continuous(name = NULL) +
  vacscale(show = FALSE) + varscale() +
  coord_cartesian(
    ylim = c(NA, 1),
    expand = FALSE
  )

veffplot <- function(
  fromdate = "2021-01-01",
  outcomes = c("c", "d"),
  flab = FALSE, xlab=TRUE
) ggplot(var.impact[date >= fromdate][var %in% outcomes]) + aes(
  date, ceff,
  linetype = c("none", "present")[mutation + 1],
  color = c("none", "present")[vac + 1],
  group = interaction(realization, vac)
) +
  facet_grid(var ~ ., scales = "free", labeller = labeller(
    var = c(c="Cumul. eff. (symptomatic)", d="Cumul. eff. (deaths)")
  ), switch = "y") +
  datebg(ylims = c(-Inf,Inf)) +
  geom_line(linetype="solid", alpha = 0.025) +
  geom_line(
    aes(group = NULL),
    data = function(dt) dt[,
                           .(ceff=median(ceff)),
                           by=.(date, mutation, var, vac)
    ]
  ) +
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
  outcomes = c("c", "d"),
  refpop = 375475/1e4,
  minrate = NA,
  flab = TRUE, xlab=FALSE
) ggplot(mlt[date >= fromdate][var %in% outcomes]) + aes(
  date, value/refpop,
  color = c("none", "present")[vac + 1],
  linetype = c("none", "present")[mutation + 1],
  group = interaction(realization, mutation, vac)
) +
  facet_grid(. ~ outcome, scales = "free_y") +
  datebg(ylims = c(0, Inf)) +
  geom_line(linetype="solid", alpha = 0.025) +
  geom_line(
    aes(group = NULL),
    data = function(dt) dt[,
      .(value=median(value)),
      by=.(date, vac, mutation, outcome)]
  ) +
  #scale_y_continuous(
  scale_y_log10(
    name = NULL,
    labels = scales::label_number_si(accuracy = 0.1)
  ) +
  vacscale() +
  varscale() +
  coord_cartesian(
    ylim = c(minrate, NA),
    expand = FALSE
  ) + (if (!flab) theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank()
  )) + (if (!xlab) theme(
    axis.text.x = element_blank()
  ))

avertplot <- function() 

outplot <- (
  epiplot() / effplot() #/ veffplot()
) + plot_layout(guides = "collect") & theme_minimal(
  base_size = 24 # base FONT size in pts
) & theme(
  legend.position = "bottom", strip.placement = "outside",
  panel.grid.major.x = element_blank()
)

ggsave(tail(.args, 1), outplot, width = 11, height = 11, units = "in", dpi = 300, bg = "white")
