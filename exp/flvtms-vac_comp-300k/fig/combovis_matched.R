
.pkgs <- c("data.table", "ggplot2", "ggrepel", "patchwork")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", c("digest.rds", "digest-key.rds", "digest-ref.rds", "vis_support.rda")),
  file.path("fig", "combo-low-matched.png")
) else commandArgs(trailingOnly = TRUE)

#' comes key'd
eff.dt <- readRDS(.args[1])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, value, averted, c.effectiveness
)][outcome %in% c("inf", "sev", "deaths")]

scn.dt <- readRDS(.args[2])
ref.dt <- readRDS(.args[3])[
  date > "2020-12-01"
][outcome %in% c("inf", "sev", "deaths")][,
  .SD, .SDcols = -c("c.value")
][, c("action", "active", "active_full", "action_full") := "none"]

load(.args[4])

mtscns <- scn.dt[
  (stockpile == "none") | # nothing & q-only
  (stockpile == "low" & active == "passive" & action == "vonly") | # as implemented
  (stockpile == "low-matched" & active == "passive" & action == "vandq") | # extra vaccine
  (stockpile == "low-matched" & active == "risk" & action == "vandq") | # risk-targetted
  (stockpile == "low" & active == "ring" & action == "vandq") | # ring targetted
  (stockpile == "covax" & action == "vandq") # LMIC, risk and ring
]

plt.dt <- setkeyv(
  {
    ret <- eff.dt[mtscns, on=.(scenario), nomatch = 0]
    ret[
      stockpile == "low-matched" & active == "passive", active := "passive+"
    ]
    ret[, stockpile := gsub("-matched", "", stockpile) ]
    ret
  },
  union(key(eff.dt), colnames(mtscns))
)

rm(eff.dt)
gc()

plt.dt[,
  action_full := fifelse(
    action == "none", as.character(action), fifelse(
    stockpile == "covax", paste(action, active, "only", sep="_"), fifelse(
    active == "none", as.character(action), paste(action, active, sep="_")
  )))
][,
  active_full := action_full
]

setkeyv(plt.dt, union(key(plt.dt), c("action_full", "active_full")))

aesbase <- aes(
  linetype = action_full, color = active_full,
  group = interaction(scenario, realization)
)

p.gen <- function(
  dt, yvar, aesbase, scale_y,
  ymax = NA, geom = geom_river
) {
  aesbase$y <- str2lang(yvar)
  aesbase$x <- quo(date)
  loc.dt <- dt[, c(key(dt), yvar), with = FALSE]
  return(ggplot(loc.dt) + aes_string(
    y = yvar
  ) +
    geom_month_background(
      loc.dt,
      font.size = 3,
      by = c(row = "outcome"), value.col = yvar, ymax = ymax
    ) +
    geom(
      mapping = aesbase
    ) +
    facet_typical(cols = NULL) +
    scale_y() +
    scale_x_null() +
    scale_color_fullscenario() +
    scale_linetype_fullscenario() +
    scale_alpha_measure(guide = "none") +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = 8, margin = margin(b = 0.1, unit = "line")), # facet
      axis.text = element_text(size = 8), # ticks
      axis.title = element_text(size = 8, margin = margin(b = 0.1, unit = "line")) # axis label
    ) + guides(
      color = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size = 2)),
      linetype = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size = 2))
    ))
}

aeseff <- aesbase
aeseff$label <- quo(sprintf('%.2f', c.effectiveness))

p.eff <- p.gen(
  plt.dt[scenario != 1],
  "c.effectiveness",
  aeseff,
  scale_y_effectiveness,
  ymax = 1
) + coord_cartesian(ylim = c(0, 1)) + theme(legend.position = "none")

p.ave <- p.gen(
  plt.dt[scenario != 1], "averted", aesbase, scale_y_averted
) + theme(legend.position = "none")

p.inc = p.gen(
  setkeyv(rbind(plt.dt[scenario != 1], ref.dt, fill = TRUE), key(plt.dt)),
  "value", aesbase, scale_y_incidence
) + theme(
  legend.position = "bottom",
  legend.margin = margin(),
  legend.spacing = unit(0.5, "line")
)

p.fin <- p.inc + p.ave + p.eff + guide_area() + plot_layout(
  design = "
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  DDD
", guides = 'collect')

ggsave(gsub("-(\\w+)\\.","-noq.",tail(.args, 1)), p.fin, height = 4*1.5, width = 3.25*4, bg = "white")

quit()

isq <- expression((action %like% "q") | (action == "none"))

p.eff <- p.gen(
  plt.dt[scenario != 1][eval(isq)], "c.effectiveness", aesbase, "Test", scale_y_effectiveness,
  ymax = 1, show.end = TRUE
) + coord_cartesian(ylim = c(0, 1)) + theme(legend.position = "none")

p.ave <- p.gen(
  plt.dt[scenario != 1][eval(isq)], "averted", aesbase, "Test", scale_y_averted
) + theme(legend.position = "none")

p.inc = p.gen(
  setkeyv(rbind(plt.dt[eval(isq)][scenario != 1], ref.dt, fill = TRUE), key(plt.dt)),
  "value", aesbase, "Test", scale_y_incidence
) + theme(
  legend.position = "bottom",
  legend.margin = margin(),
  legend.spacing = unit(0.5, "line")
)

p.fin <- p.inc + p.ave + p.eff + guide_area() + plot_layout(
  design = "
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  ABC
  DDD
", guides = 'collect')

ggsave(gsub("-(\\w+)\\.","-isq.",tail(.args, 1)), p.fin, height = 4*1.5, width = 3.25*4, bg = "white")
