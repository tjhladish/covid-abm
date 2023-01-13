
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest-ref.rds", "digest-key.rds", "vocwindows.rds", "digest.rds")),
  file.path("fig", "output", "alt_inc_non.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

intfilter <- expression(realization >= 0)

#' comes key'd - all the intervention results
inc.dt <- readRDS(.args[2])[
  eval(datefilter) & eval(intfilter) & eval(outfilter)
][, .(
  scenario, realization, date, outcome, value
)]

scn.dt <- readRDS(.args[3])
takeover.wins <- readRDS(.args[4])
# reconstructing reference scenarios
refscn.dt <- scn.dt[(quar == FALSE & act_vac == "none") & (!inf_con | !pas_vac)]

int.dt <- readRDS(.args[5])[
  eval(datefilter) & eval(intfilter) & eval(outfilter)
][scenario %in% refscn.dt$scenario, .(
  scenario, realization, date, outcome, value
)]

plt.dt <- setkeyv(
  rbind(inc.dt, int.dt)[refscn.dt, on=.(scenario)],
  union(key(inc.dt), colnames(scn.dt))
)

rm(inc.dt)
gc()

plt.qs <- plt.prep(plt.dt, j = .(value))

mm.ref <- plt.qs[,.(ymin = min(q90l), ymax = max(q90h)),by=.(outcome)]
mm.ref[, ymin := ymin - .15*(ymax-ymin) ]
mm.ref[, yspan := ymax - ymin ]
tw <- takeover.wins[q == 0.5][CJ(measure, outcome = mm.ref$outcome), on=.(measure)][mm.ref, on=.(outcome)]
tw[, end := pmin(end, plt.qs[, max(date)])]
tw[, mid := start + (end-start)/2 ]

allplot <- function (data.qs, yl, withRef = FALSE, col.breaks = if (withRef) c("risk",
                                                                    "age", "ring") else c("none", "risk", "age", "ring"), withBands = NULL,
          ins = list())
{
  res <- ggplot(data.qs) + aes(
    x = date, color = pas_alloc, linetype = factor(c("nonpi",
                                                                              "wquar")[quar + 1])) + facet_nested(rows = vars(outcome),
                                                                                                                   switch = "y", scales = "free_y",
                                                                                                                  labeller = labeller(outcome = c(inf = "Infection", sev = "Severe Disease",
                                                                                                                                                  deaths = "Deaths", vaccine = "Per 10K,\nDoses Administered"))) +
    geom_month_background(data.qs, by = c(row = "outcome",
                                          col = "talloc"), font.size = 3, value.col = "qmed",
                          max.col = "q90h", min.col = "q90l") + ins
  if (withRef) {
    res <- res + geom_hline(aes(yintercept = 0, color = "none",
                                linetype = "nonpi"), show.legend = FALSE, data = function(dt) dt[,
                                                                                                 .SD[1], by = .(outcome, talloc)])
  }
  res <- res + geom_ribbon(aes(ymin = q90l, ymax = q90h, fill = pas_alloc,
                               color = NULL), alpha = 0.15) + geom_line(aes(y = qmed)) +
    scale_color_discrete(name = "Supply", aesthetics = c("color", "fill")) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(name = yl, expand = c(0, 0), labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_x_null() + scale_linetype_quar(guide = "none") +
    theme_minimal() +
    theme(legend.position = "bottom", strip.placement = "outside",
          legend.direction = "horizontal")
  return(res)
}

p <- allplot(
  plt.qs, yl = "Per 10k, Incidence of ...",
  ins = voc.band(tw)
)

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
