
#' DANGER
#' because this ends with saving everything, have to clean env first when
#' used interactively.
#' TODO: do via creating and saving an environment, rather than everything
rm(list = ls())

.pkgs <- c("data.table", "ggplot2", "ggrepel")

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
  name = "Cumulative Effectiveness Against ...",
  breaks = seq(0, 1, by=0.25)
)

scale_y_doses <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Per 10k, Daily Doses Administered"
)

scale_y_cdoses <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Per 10k, Cumulative Doses Administered"
)

scale_y_incidence <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Per 10k, Incidence of ..."
)

scale_y_averted <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Per 10k, Averted Incidence of ..."
)

scale_x_null <- gg_scale_wrapper(
  scale_x_continuous,
  name = NULL, breaks = NULL, minor_breaks = NULL,
  expand = expansion(
    mult = c(0, 0), add = c(0, 0)
  )
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
  breaks = c("none", "passive", "ring", "risk"),
  labels = c(none = "None", passive = "Mass Only", ring = "Mass+Ring", risk = "Mass+High Risk"),
  values = c(none = "black", passive = "#fb6502",  ring = "#00529b", risk = "#006b35"),
  drop = TRUE, limits = force
)

fullscnbreaks <- c(
  "none", "qonly",
  "vonly_passive", "vandq_passive",
  "vonly_passive+", "vandq_passive+",
  "vonly_risk", "vandq_risk",
  "vonly_ring", "vandq_ring"
)

fullscnlabels <- c(
  vandq_passive ="Mass Vac. & Con. Tracing => Quar.",
  vonly_passive = "Mass Vaccination",
  "vandq_passive+" ="Expanded Mass Vac. & Con. Tracing => Quar.",
  "vonly_passive+" = "Expanded Mass Vac.",
  vandq_risk = "Mass + High Risk Vac. & Con. Tracing => Quar.",
  vonly_risk = "Mass + High Risk Vac.",
  vandq_ring="Mass Vac. & Con. Tracing => Ring Vac. + Quar.",
  vonly_ring = "Mass Vac. & Con. Tracing => Ring Vac.",
  qonly = "Contact Tracing => Quarantine",
  none = "None"
)

fullscnlines <- c(
  vandq_ring="solid",
  vandq_passive = "solid",
  "vandq_passive+" = "solid",
  vandq_risk = "solid",
  vonly_ring = "dashed",
  vonly_risk = "dashed",
  vonly_passive = "dashed",
  "vonly_passive+" = "dashed",
  qonly = "dotted", none = "solid"
)

fullscncolors <- c(
  vandq_ring = "#00529b",
  vonly_ring = "#00529b",
  vandq_passive="#fb6502",
  vonly_passive = "#fb6502",
  "vandq_passive+" = "#fb6502",
  "vonly_passive+" = "#fb6502",
  vandq_risk = "#006b35",
  vandq_risk = "#006b35",
  qonly = "black",
  none = "black"
)

scale_linetype_fullscenario <- gg_scale_wrapper(
  scale_linetype_manual,
  name = "Intervention",
  breaks = fullscnbreaks,
  labels = fullscnlabels,
  values = fullscnlines,
  drop = TRUE, limits = force
#  , guide = guide_legend(override.aes = list(color = "grey"))
)

scale_color_fullscenario <- gg_scale_wrapper(
  scale_color_manual,
  name = "Intervention",
  breaks = fullscnbreaks,
  labels = fullscnlabels,
  values = fullscncolors,
  drop = TRUE, limits = force
  #  , guide = guide_legend(override.aes = list(color = "grey"))
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
    value.col = "value",
    by = setdiff(key(dt), c("realization", "date")),
    #' TODO by doesn't currently handle nested facets
    date.col = "date",
    mcycle = c("off", "on"),
    ymax = NA
) {
  dt <- as.data.table(dt)
  if ("col" %in% names(by)) {
    #' need to consider ymax across rows if by["col"] exists
    ymref <- dt[, c(max(get(value.col), ymax, na.rm = TRUE)), by = eval(unname(by["row"]))]
    dt[ymref, ym := V1, on=unname(by["row"])]
    dt <- dt[!is.na(get(by["col"]))]
  } else {
    dt[, ym := max(get(value.col), ymax, na.rm = TRUE), by = eval(unname(by)) ]
  }
  dt[, {
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
      mtype = rep(mcycle, length.out = length(starts)),
      mon = month(starts),
      yr = year(starts),
      yshow = c(TRUE, month(starts[-1])==1),
      ymax = ym[1]
    )
  }, by = eval(unname(by))]
}

m.abb <- gsub("^(.).+$","\\1", month.abb)

#' then:
#'  - paint month-wide bars, starting with first full enclosed month
#'  - label months from lower bound to upper bound month at midpoints
#'  - label first fully enclosed month
geom_month_background <- function(
    data,
    col.cycle = c(off = NA, on = alpha("lightgrey", 0.75)),
    m.labels = m.abb,
    font.size = 5,
    ylog = FALSE,
    datafn = tsref.dt,
    ...
) {
  text.col <- alpha(col.cycle, 1)
  names(text.col) <- c(tail(names(col.cycle), -1), names(col.cycle)[1])
  text.col[is.na(text.col)] <- "white"
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]
  dt[, col := text.col[mtype] ]
  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (ylog) 0 else -Inf, ymax = Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    ),
    geom_text(
      mapping = aes(x = mid, y = ymax*.9, label = m.abb[mon]),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, color = dt$col,
      size = font.size, vjust = "bottom"
    ),
    geom_text(
      mapping = aes(x = mid, y = ymax*.85, label = yr),
      data = dt[yshow == TRUE], angle = 90,
      inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
      size = font.size, hjust = "right"
    )
    # ,
    # geom_text(
    #   mapping = aes(x = mid, y = 0.25, label = yr),
    #   data = dt[yshow == TRUE],
    #   inherit.aes = FALSE, show.legend = FALSE, color = rep(dt[yshow == TRUE]$col, facet.mul)
    # )
  )
}

med.dt <- function(
  dt,
  yvar
) {
  bynames <- setdiff(
    colnames(dt),
    c("realization", yvar, "value", "averted", "c.averted", "c.effectiveness", "c0.effectiveness")
  )
  setnames(dt[,.(median(get(yvar))), keyby = bynames ], "V1", yvar)
}

geom_spaghetti <- function(
  mapping, data,
  alpha = 0.01,
  show.end = FALSE,
  spag.size = 0.1,
  max.lines = if (interactive()) 10 else 150,
  sample.var = "realization",
  ...
) {
  aesmany <- mapping
  aesmany$linetype <- NULL
  aesone <- mapping
  aesone$group <- quote(scenario)
  mdt <- med.dt(data, rlang::as_name(aesmany$y))
  ret <- list(
    geom_line(aesmany, size = spag.size, data = data[get(sample.var) < max.lines], alpha = alpha, ...),
    geom_line(aesone, data = mdt, ...)
  )
  if (show.end) {
    aesend <- aesmany
    aesend$group <- aesone$group
    aesend$label <- str2lang(paste0("sprintf('%.2f',", rlang::as_name(aesend$y),")"))
    ret[[length(ret)+1]] <- geom_text_repel(
      aesend, data = mdt[date == max(date)],
      hjust = "left", direction = "y", nudge_x = 7,
      size = 2, segment.size = 0.1, min.segment.length = 0,
      show.legend = FALSE
    )
  }
  ret
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
    stockpile = c(
      low = "(low) XXk per 10k-day", high = "(high) YYk per 10k-day",
      "low-matched" = "XXk per 10k-day\nExpanded Mass Vac.",
      "high-matched" = "YYk per 10k-day\nExpanded Mass Vac.")
  )
)

save(list=ls(), file = tail(.args, 1))
