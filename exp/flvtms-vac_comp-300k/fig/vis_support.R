
#' DANGER
#' because this ends with saving everything, have to clean env first when
#' used interactively.
#' TODO: do via creating and saving an environment, rather than everything
rm(list = ls())

.pkgs <- c("data.table", "ggplot2")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args <- if (interactive()) {
  file.path("fig", "vis_support.rda")
} else commandArgs(trailingOnly = TRUE)

gg_scale_wrapper <- function(
    scale_fun,
    ...
) {
  stopifnot(!missing(scale_fun))
  defs <- list(...)
  if (!length(defs)) warning(
    "provided no default arguments; consider using scale_fun directly."
  )

  return(function(...) {
    # this different ... is how you get a function back that let's you
    # override defaults, set other arguments to scale_... functions
    .ellipsis <- list(...)
    .args <- defs
    .args[names(.ellipsis)] <- .ellipsis
    do.call(scale_fun, .args)
  })
}

scale_y_effectiveness <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Cumulative Effectiveness Against...",
  breaks = seq(0, 1, by=0.25)
)

scale_y_incidence <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Per 10k (log2 scale), Incidence of ...",
  trans = "log2"
)

scale_y_averted <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Per 10k, Averted Incidence of ..."
)

scale_x_null <- gg_scale_wrapper(
  scale_x_continuous,
  name = NULL, breaks = NULL, minor_breaks = NULL
)

scale_linetype_scenario <- gg_scale_wrapper(
  scale_linetype_manual,
  name = "Intervention",
  breaks = c("none", "qonly", "vonly", "vandq"),
  labels = c(vandq="Vac.+Quar.", vonly = "Vac. Only", qonly = "Quar. Only", none = "None"),
  values = c(vandq="solid", vonly = "dashed", qonly = "dotted", none = "solid"),
  drop = TRUE, limits = force,
  guide = guide_legend(override.aes = list(color = "grey"))
)

scale_color_scenario <- gg_scale_wrapper(
  scale_color_manual,
  name = "Vaccine Distr.",
  breaks = c("none", "passive", "ring"),
  labels = c(none = "None", passive = "Mass Only", ring = "Mass+Ring"),
  values = c(none = "black", passive = "blue",  ring = "green"),
  drop = TRUE, limits = force
)

#' @title extract month bounds
#'
#' @description takes time series data and extracts the bounding month windows,
#' i.e. 1st of month to last of month, as an ordered series of rows in
#' \code{\link[data.table]{data.table}}.
#'
#' @param dt an object cast-able via \code{\link[data.table]{as.data.table}}
#' @param date.col character; the column corresponding to the date information;
#'   column should be cast-able via \code{\link[base]{as.Date}}
#' @param mcycle the month cycle labeling - pattern is repeated for length of
#'   resulting bounds
#'
#' @details This method takes a `data.table`-compatible structure with a `Date`-
#' compatible column, and returns a `data.table` with rows indicating the
#' bounding month start & end dates, plus cyclical labels for those rows,
#' and additional pre-calculated
#'
#' @return a `data.table`, columns `start` and `end` (both `Date`),
#' `mtype` (`factor`), `mon` (`integer`, 1-12),
#' `yshow` (`logical`, `TRUE` for first entry & Januarys)
#'
#' examine a time series of data, determine
#'  - the first full enclosed month
#'  - the lower bound month (if different)
#'  - the upper bound month
tsref.dt <- function(
    dt,
    date.col = "date",
    mcycle = c("off", "on")
) {
  as.data.table(dt)[, {
    # assert: guarantees YYYY-MMM-01 to YYYY-MMM-01
    date.range <- trunc(range(as.Date(get(date.col))), "months")
    date.min <- date.range[1]
    # goes into next month, truncates back to 1st, then ticks back to last of original month
    date.max <- trunc(date.range[2] + 35, "months")-1
    # extract the inclusive month window start/ends
    starts <- seq(date.min, date.max, by="months")
    ends <- c(starts[-1] - 1, date.max)
    .(
      start = starts,
      end = ends,
      mid = as.numeric(ends - starts)/2 + starts,
      mtype = mcycle,
      mon = month(starts),
      yr = year(starts),
      yshow = c(TRUE, month(starts[-1])==1)
    )
  }]
}

m.abb <- gsub("^(.).+$","\\1", month.abb)

#' then:
#'  - paint month-wide bars, starting with first full enclosed month
#'  - label months from lower bound to upper bound month at midpoints
#'  - label first fully enclosed month
geom_month_background <- function(
    data,
    facet.mul,
    col.cycle = c(off = NA, on = alpha("lightgrey", 0.75)),
    m.labels = m.abb,
    font.size = 5,
    ylog = FALSE
) {
  text.col <- alpha(col.cycle, 1)
  names(text.col) <- c(tail(names(col.cycle), -1), names(col.cycle)[1])
  text.col[is.na(text.col)] <- "white"
  dt <- tsref.dt(data)
  dt[, fill := col.cycle[mtype] ]
  dt[, col := text.col[mtype] ]
  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (ylog) 0 else -Inf, ymax = Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = rep(dt$fill, facet.mul)
    ),
    geom_text(
      mapping = aes(x = mid, y = 0.5, label = m.abb[mon]),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, color = rep(dt$col, facet.mul),
      size = font.size, vjust = "bottom"
    ),
    geom_text(
      mapping = aes(x = mid, y = 0.5-0.05, label = yr),
      data = dt[yshow == TRUE], angle = -90,
      inherit.aes = FALSE, show.legend = FALSE, color = rep(dt[yshow == TRUE]$col, facet.mul),
      size = font.size, hjust = "left"
    )
    # ,
    # geom_text(
    #   mapping = aes(x = mid, y = 0.25, label = yr),
    #   data = dt[yshow == TRUE],
    #   inherit.aes = FALSE, show.legend = FALSE, color = rep(dt[yshow == TRUE]$col, facet.mul)
    # )
  )
}

med.dt <- function(dt) {
  bynames <- setdiff(
    colnames(dt),
    c("realization", "value", "averted", "c.averted", "c.effectiveness")
  )
  tar <- intersect(colnames(dt),c("value", "averted", "c.averted", "c.effectiveness"))[1]
  setnames(dt[,.(median(get(tar))), keyby = bynames ], "V1", tar)
}

geom_spaghetti <- function(
  mapping, data,
  alpha = 0.02
) {
  aesmany <- mapping
  aesmany$linetype <- NULL
  aesone <- mapping
  aesone$group <- quote(scenario)
  list(
    geom_line(aesmany, size = 0.05, data = data, alpha = alpha),
    geom_line(aesone, data = med.dt(data))
  )
}

gg_facet_wrapper <- function(
    facet_fun,
    ...
) {
  stopifnot(!missing(facet_fun))
  defs <- list(...)
  if (!length(defs)) warning(
    "provided no default arguments; consider using facet_fun directly."
  )

  return(function(...) {
    # this different ... is how you get a function back that let's you
    # override defaults, set other arguments to scale_... functions
    .ellipsis <- list(...)
    .args <- defs
    .args[names(.ellipsis)] <- .ellipsis
    do.call(facet_fun, .args)
  })
}

facet_typical <- gg_facet_wrapper(
  facet_grid,
  rows = vars(outcome),
  cols = vars(stockpile),
  scales = "free_y",
  drop = TRUE,
  switch = "y",
  labeller = labeller(
    outcome = c(inf = "Infection", symp = "Symptoms", sev = "Severe Disease", crit = "Critical Disease", deaths = "Death"),
    stockpile = c(low = "(low) XXk per 10k-day", high = "(high) YYk per 10k-day")
  )
)

save(list=ls(), file = tail(.args, 1))
