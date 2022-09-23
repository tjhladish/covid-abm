
.pkgs <- c(
  "data.table", "ggplot2", "scales",
  "ggh4x", "cabputils", "patchwork"
)

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- commandArgs(
  trailingOnly = TRUE,
  args = c(
    file.path("fig", "vis_support.rda"),
    file.path("fig", "process", c("digest.rds", "digest-key.rds")),
    file.path("fig", "combo_alt.png")
  )
)

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  eval(datefilter)
][, .(
  scenario, realization, date, outcome, value, averted, c.effectiveness
)]

eff.dt[, c.averted := cumsum(averted), by=setdiff(key(eff.dt), "date")]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

#' TODO key order?
plt.dt <- eff.dt[scn.dt, on=.(scenario)] |> melt.data.table(
  id.vars = union(key(eff.dt), colnames(scn.dt)),
  variable.name = "measure"
) |> setkeyv(union(key(eff.dt), c(colnames(scn.dt), "measure"))) |>
  quantile.data.table(
    samplecol = "realization",
    probs = qprobs(c(`50`=.50, `90`=.90), TRUE, FALSE)
  )

rm(eff.dt)
gc()

scale_color_state <- gg_scale_wrapper(
  scale_color_manual,
  name = NULL,
  breaks = c("MS", "FL", "VT"),
  values = c("#F8766D", "#00BA38", "#619CFF"),
  aesthetics = c("color", "fill")
)

p.inc <- ggplot(plt.dt[measure != "averted"]) + aes(
  x = date
) +
  facet_grid2(
    outcome ~ measure, scales = "free_y", independent = "y",
    labeller = labeller(
      outcome = c(inf = "Infection", sev = "Severe Disease", deaths = "Death"),
      measure = c(
        "value" = "Per 10k,\nIncidence of ...",
        "c.averted" = "Per 10k,\nCumulative Averted ...",
        "c.effectiveness" = "Cumulative Effectiveness\nAgainst..."
      )
    ), switch = "y", drop = TRUE
  ) +
  geom_month_background(
    plt.dt[measure != "averted"],
    by = c("measure", "outcome"),
    font.size = 3, value.col = "q90h"
  ) +
  geom_ribbon(
    aes(ymin=q90l, ymax=q90h, fill = state),
    alpha = 0.25
  ) +
  geom_line(
    aes(y=qmed, color=state)
  ) +
  scale_x_null() +
  scale_color_state() +
  scale_y_continuous(name = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

# p.eff <- ggplot(plt.dt[measure == "c.effectiveness"]) + aes(
#   x = date
# ) +
#   facet_nested(
#     outcome ~ measure,
#     labeller = labeller(
#       outcome = c(inf = "Infection", sev = "Severe Disease", deaths = "Death"),
#       measure = c("c.effectiveness" = "Cumulative Effectiveness\nAgainst...")
#     ), drop = TRUE
#   ) +
#   geom_month_background(
#     plt.dt[measure == "c.effectiveness"],
#     by = c("measure", "outcome"),
#     font.size = 3, value.col = "q90h"
#   ) +
#   geom_ribbon(
#     aes(ymin=q90l, ymax=q90h, fill = state),
#     alpha = 0.25
#   ) +
#   geom_line(
#     aes(y=qmed, color=state)
#   ) +
#   scale_x_null() +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_color_state() +
#   scale_y_continuous(name = NULL, position = "right") +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom", strip.placement = "outside"
#   )
#
# p <- p.inc + p.eff

ggsave(tail(.args, 1), p.inc, height = 6, width = 10, bg = "white")
