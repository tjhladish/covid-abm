
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

scale_y_fraction <- gg_scale_wrapper(
  scale_y_continuous,
  name = "Fraction of Population",
  breaks = (0:4)/4,
  limits = c(0, 1),
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
  drop = TRUE, limits = force,
  aesthetics = c("color", "fill")
)

fullscnbreaks <- c(
  "none", "qonly",
  "vonly_passive", "vandq_passive",
  "vonly_passive+", "vandq_passive+",
  "vonly_risk", "vandq_risk",
  "vonly_ring", "vandq_ring",
  "vonly_risk_only", "vandq_risk_only",
  "vonly_ring_only", "vandq_ring_only"
)

fullscnlabels <- c(
  vandq_passive ="Vac. & Con. Tracing => Quar.",
  vonly_passive = "Vaccination",
  "vandq_passive+" ="Expanded Vac. & Con. Tracing => Quar.",
  "vonly_passive+" = "Expanded Vac.",
  vandq_risk = "Vac + High Risk & Con. Tracing => Quar.",
  vonly_risk = "Vac + High Risk",
  vandq_ring="Vac. & Con. Tracing => Ring Vac. + Quar.",
  vonly_ring = "Vac. & Con. Tracing => Ring Vac.",
  qonly = "Contact Tracing => Quarantine",
  none = "None",
  "vandq_risk_only" = "High Risk Vac. & Con. Tracing => Quar.",
  "vonly_risk_only" = "High Risk Vac.",
  "vandq_ring_only" = "Con. Tracing => Ring Vac. + Quar.",
  "vonly_ring_only" = "Con. Tracing => Ring Vac."
)

fullscnlines <- c(
  vandq_ring="solid",
  vandq_passive = "solid",
  "vandq_passive+" = "solid",
  vandq_risk = "solid",
  vonly_ring = "dashed",
  vonly_risk = "dashed",

  "vandq_ring_only" = "solid",
  "vandq_risk_only" = "solid",
  "vonly_ring_only" = "dashed",
  "vonly_risk_only" = "dashed",

  vonly_passive = "dashed",
  "vonly_passive+" = "dashed",
  qonly = "dotted", none = "solid"
)

fullscncolors <- c(
  vandq_ring = "#00529b",
  vonly_ring = "#00529b",
  "vandq_ring_only" = "firebrick",
  "vonly_ring_only" = "firebrick",
  "vandq_risk_only" = "dodgerblue",
  "vonly_risk_only" = "dodgerblue",

  vandq_passive="#fb6502",
  vonly_passive = "#fb6502",
  "vandq_passive+" = "#fb6502",
  "vonly_passive+" = "#fb6502",
  vandq_risk = "#006b35",
  vonly_risk = "#006b35",
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
  , guide = guide_legend(override.aes = list(fill = NA))
)

scale_color_fullscenario <- gg_scale_wrapper(
  scale_color_manual,
  name = "Intervention",
  breaks = fullscnbreaks,
  labels = fullscnlabels,
  values = fullscncolors,
  drop = TRUE, limits = force,
  aesthetics = c("color", "fill")
  , guide = guide_legend(override.aes = list(fill = NA))
)

scale_color_measure <- gg_scale_wrapper(
  scale_color_manual,
  name = NULL,
  values = c(
    case = "royalblue3", death = "green4", infection = "orangered",
    hospPrev = "orange", hospInc = "orange", vaxHosp = 'dodgerblue4',
    brkthru = "tan4"
  ),
  guide = "none"
)

scale_color_inputs <- gg_scale_wrapper(
  scale_color_manual,
  name = NULL,
  values = c(
    socialdist = "darkorange3", seasonality = "purple",
    vocprev1 = 'royalblue3', vocprev2 = 'turquoise4', vocprev3 = 'darkorchid3',
    coverage = 'purple'
  ),
  guide = "none",
  aesthetics = c("color", "fill")
)

measlbls <- c(
  observed = "Observed", sample = "Simulated Replicates",
  central = "Median Simulation",
  quantile = "90% Prediction Interval"
)

scale_alpha_measure <- gg_scale_wrapper(
  scale_alpha_manual,
  name = NULL,
  breaks = c("observed", "sample", "central", "quantile"),
  labels = measlbls,
  values = c(observed = 0.6, sample = 0.05, quantile = 0.5, central = 1),
  guide = guide_legend(
    override.aes = list(
      linetype = c("blank", "solid", "solid"),
      size = c(10, 1, 3)
    )
  )
)

scale_shape_measure <- gg_scale_wrapper(
  scale_shape_manual,
  name = NULL,
  breaks = c("observed", "sample", "central"),
  labels = measlbls,
  values = c(observed = 20, sample = NA, central = NA),
  guide = guide_legend(
    override.aes = list(linetype = c("blank", "solid", "solid"))
  )
)

geom_observation <- gg_scale_wrapper(
  geom_point,
  mapping = aes(alpha = "observed", shape = "observed"),
  data = function(dt) subset(dt, is.na(realization)),
  size = 2.5,
  stroke = 0
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
    ymax = NA, ymin = NA
) {
  dt <- as.data.table(dt)
  maxer <- function(x) min(max(x, ymax, na.rm = TRUE), ymax, na.rm = TRUE)
  miner <- function(x) max(min(x, ymin, na.rm = TRUE), ymin, na.rm = TRUE)
  if(!is.null(by)) { if ("col" %in% names(by)) {
    #' need to consider ymax across rows if by["col"] exists
    ymref <- dt[, .(mx = maxer(get(value.col)), mn = miner(get(value.col))), by = eval(unname(by["row"]))]
    dt[ymref, c("ym", "yn") := .(mx, mn), on=unname(by["row"])]
    dt <- dt[!is.na(get(by["col"]))]
  } else {
    dt[, c("ym","yn") := .(
      maxer(get(value.col)),
      miner(get(value.col))
    ), by = eval(unname(by)) ]
  } } else {
    dt[, c("ym","yn") := .(
      maxer(get(value.col)),
      miner(get(value.col))
    ) ]
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
      ymax = ym[1], ymin = yn[1]
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
    font.size = 8, font.face = "bold",
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
  xform <- if (ylog) {
    function(hi, lo, dropto) lo*(hi/lo)^dropto
  } else {
    function(hi, lo, dropto) lo+(hi-lo)*dropto
  }
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (ylog) 0 else -Inf, ymax = Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    ),
    geom_text(
      mapping = aes(x = mid, y = xform(ymax, ymin, .87), label = m.abb[mon]),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, color = dt$col,
      size = font.size, vjust = "bottom", fontface = font.face
    ),
    geom_text(
      mapping = aes(x = mid, y = xform(ymax, ymin, .83), label = yr),
      data = dt[yshow == TRUE], angle = 90,
      inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
      size = font.size, hjust = "right", fontface = font.face
    )
    # ,
    # geom_text(
    #   mapping = aes(x = mid, y = 0.25, label = yr),
    #   data = dt[yshow == TRUE],
    #   inherit.aes = FALSE, show.legend = FALSE, color = rep(dt[yshow == TRUE]$col, facet.mul)
    # )
  )
}


#' @title Quantile [data.table] Groups
#'
#' @description Applies [quantile()] to [data.table] groups,
#' across sample column, for target column(s). Returns results by
#' group and (when more than one target column) measure.
#'
#' @param dt a [data.table]
#'
#' @param col the column(s) to quantile
#'
#' @param samplecol the sample-defining column (i.e., what gets
#' quantiled over)
#'
#' @param keyby the grouping definition
#'
#' @param ... arguments to [stats::quantile()]. If `probs` provided
#' with names, those names will be used for the resulting quantile columns
#'
#' @return a [data.table], cols for:
#'  * the keys
#'  * optionally, `measure` if `length(col) > 1`
#'  * a column for each quantile returned; if `probs` in `...` and named,
#'  those are the column names
quantile.data.table <- function(
  dt, col = "value",
  samplecol,
  keyby = setdiff(key(dt), samplecol),
  ...
) {

  if (!all(col %in% colnames(dt))) {
    stop("`colnames(dt)` and `col` mismatch")
  }

  if (length(col) > 1) {
    idcol <- "measure"
    names(col) <- col
  } else {
    idcol <- NULL
  }

  dots <- list(...)
  curryq <- if (!is.null(dots$probs) && !is.null(names(dots$probs))) {
    function(x) quantile |> do.call(c(list(x), dots)) |> setNames(names(dots$probs))
  } else {
    function(x) quantile |> do.call(c(list(x), dots))
  }

  return(
    col |> lapply(\(tarcol) {
      dt[,{
        get(tarcol) |> curryq() |> as.list()
      }, keyby = keyby]
    }) |> rbindlist(idcol = idcol) |>
    setkeyv(cols = c(idcol, keyby))
  )

}

geom_spaghetti <- function(
  mapping,
  data,
  spaghetti = {
    grouptree <- as.character(rlang::get_expr(mapping$group))
    res <- list(sample = tail(grouptree, 1))
    if (grouptree[1] == "interaction") {
      if (length(grouptree) == 3) { # group = interaction(something, samplevar)
        res$subgroup <- rlang::parse_expr(grouptree[2])
      } else { # group = interaction(...somethings..., samplevar)
        res$subgroup <- rlang::parse_expr(head(grouptree, -1))
      }
    } else { # group = samplevar
      res$subgroup <- NULL
    }
    res
  },
  datafn = \(dt) dt |> quantile(
    col = rlang::as_name(mapping$y),
    samplecol = spaghetti$sample,
    probs = c(md = 0.5)
  ) |> setnames("md", rlang::as_name(mapping$y)),
  sample.size = 0.1,
  center.size = 3*sample.size,
  max.lines = if (interactive()) 10 else 150,
  geom = geom_line,
  ...
) {
  show.end <- !is.null(mapping$label)
  # split mapping ...
  aesmany <- mapping; aesone <- mapping
  # ignore specific elements
  aesmany$linetype <- aesmany$label <- aesone$label <- NULL

  if(is.null(mapping$alpha)) {
    aesmany$alpha <- "sample"
    aesone$alpha <- "central"
  }

  if (!is.null(spaghetti$subgroup)) {
    aesone$group <- spaghetti$subgroup
  }

  if (!is.null(max.lines)) {
    datamany <- \(dt) dt[get(spaghetti$sample) < max.lines]
  } else {
    datamany <- NULL
  }
  dataone <- datafn

  if (!missing(data)) {
    dataone <- dataone(data)
    if (!is.null(datamany)) datamany <- datamany(data)
  }

  ret <- list(
    geom(aesmany, data = datamany, size = sample.size, ...),
    geom(aesone, data = dataone, size = center.size, ...)
  )

  if (show.end) {
    aesend <- mapping
    aesend$group <- aesone$group
    aesend$alpha <- aesone$alpha
    aesend$linetype <- NULL

    if (is.function(dataone)) {
      datalbl <- \(dt) {
        res <- dataone(dt)
        res[order(get(rlang::as_name(aesend$x))),
          .SD[.N], by = setdiff(key(res), rlang::as_name(aesend$x))
        ]
      }
    } else {
      datalbl <- dataone[
        order(get(rlang::as_name(aesend$x))),
        .SD[.N],
        by = setdiff(key(dataone), rlang::as_name(aesend$x))
      ]
    }

    # aesend$label <- str2lang(paste0("sprintf('%.2f',", rlang::as_name(aesend$y),")"))

    ret[[length(ret)+1]] <- geom_text_repel(
      aesend, data = datalbl,
      hjust = "left", direction = "y",
      nudge_x = 7, size = 2,
      segment.size = sample.size,
      min.segment.length = 0,
      show.legend = FALSE
    )
  }
  ret
}

geom_river <- function(
    mapping,
    data,
    spaghetti = {
      grouptree <- as.character(rlang::get_expr(mapping$group))
      res <- list(sample = tail(grouptree, 1))
      if (grouptree[1] == "interaction") {
        if (length(grouptree) == 3) { # group = interaction(something, samplevar)
          res$subgroup <- rlang::parse_expr(grouptree[2])
        } else { # group = interaction(...somethings..., samplevar)
          res$subgroup <- rlang::parse_expr(head(grouptree, -1))
        }
      } else { # group = samplevar
        res$subgroup <- NULL
      }
      res
    },
    datafn = \(dt) dt |> quantile(
      col = rlang::as_name(mapping$y),
      samplecol = spaghetti$sample,
      probs = c(lo = 0.05, md = 0.5, hi = 0.95)
    ) |> setnames("md", rlang::as_name(mapping$y)),
    sample.size = 0.1,
    center.size = 3*sample.size,
    ...
) {
  show.end <- !is.null(mapping$label)
  # split mapping ...
  aesmany <- mapping; aesone <- mapping
  # ignore specific elements
  aesmany$linetype <- aesmany$label <- aesone$label <- NULL
  aesmany$fill <- aesmany$colour
  aesmany$colour <- NULL

  # TODO generify
  aesmany$ymax <- quo(hi)
  aesmany$ymin <- quo(lo)
  aesmany$y <- NULL

  if(is.null(mapping$alpha)) {
    aesmany$alpha <- "quantile"
    aesone$alpha <- "central"
  }

  if (!is.null(spaghetti$subgroup)) {
    aesmany$group <- aesone$group <- spaghetti$subgroup
  } else {
    aesmany$group <- aesone$group <- NULL
  }

  dataone <- datafn

  if (!missing(data)) {
    dataone <- dataone(data)
  }

  ret <- list(
    geom_ribbon(aesmany, data = dataone, ...),
    geom_line(aesone, data = dataone, size = center.size, ...)
  )

  if (show.end) {
    aesend <- mapping
    aesend$group <- aesone$group
    aesend$alpha <- aesone$alpha
    aesend$linetype <- NULL

    if (is.function(dataone)) {
      datalbl <- \(dt) {
        res <- dataone(dt)
        res[order(get(rlang::as_name(aesend$x))),
            .SD[.N], by = setdiff(key(res), rlang::as_name(aesend$x))
        ]
      }
    } else {
      datalbl <- dataone[
        order(get(rlang::as_name(aesend$x))),
        .SD[.N],
        by = setdiff(key(dataone), rlang::as_name(aesend$x))
      ]
    }

    # aesend$label <- str2lang(paste0("sprintf('%.2f',", rlang::as_name(aesend$y),")"))

    ret[[length(ret)+1]] <- geom_text_repel(
      aesend, data = datalbl,
      hjust = "left", direction = "y",
      nudge_x = 7, size = 2,
      segment.size = sample.size,
      min.segment.length = 0,
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

geom_crosshair <- function(mapping, data, ...) {
  #'
  mapv <- mapping
  mapv$xmax <- NULL; mapv$xmin <- NULL
  mapv$shape <- NULL
  maph <- mapping
  maph$ymax <- NULL; maph$ymin <- NULL
  maph$shape <- NULL
  showpoint <- !is.null(mapping$shape)
  res <- list(
    geom_linerange(mapv, data = data, show.legend = FALSE),
    geom_linerange(maph, data = data, show.legend = FALSE)
  )
  if (showpoint) {
    mapp <- mapping
    mapp$xmax <- NULL
    mapp$xmin <- NULL
    mapp$ymax <- NULL
    mapp$ymin <- NULL
    mapp$alpha <- mapp$shape
    res[[length(res)+1]] <- geom_observation(mapping = mapp, data = seroprev)
  }
  return(res)
}

prepare <- function(...) setkey(melt(
  rbind(..., fill = TRUE),
  id.vars = c("realization", "date"),
  variable.name = "measure", variable.factor = FALSE
), measure, realization, date)

save(list=ls(), file = tail(.args, 1))

# TODO figure out if this could work?
# geom_month_background_alt <- function(
    #     data,
#     col.cycle = c(off = NA, on = alpha("lightgrey", 0.75)),
#     labels = m.abb,
#     font.size = 5,
#     ylog = FALSE,
#     datafn = tsref.dt,
#     ...
# ) {
#   text.col <- alpha(col.cycle, 1)
#   names(text.col) <- c(tail(names(col.cycle), -1), names(col.cycle)[1])
#   text.col[is.na(text.col)] <- "white"
#   dt <- datafn(data, ...)
#   dt[, fill := col.cycle[mtype] ]
#   dt[, col := text.col[mtype] ]
#
#   list(
#     geom_rect(
#       mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = after_scale(0), ymax = after_scale(1)),
#       data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
#     ),
#     geom_text(
#       mapping = aes(x = mid, y = after_scale(.9), label = labels[mon]),
#       data = dt, inherit.aes = FALSE, show.legend = FALSE, color = dt$col,
#       size = font.size, vjust = "bottom"
#     ),
#     geom_text(
#       mapping = aes(x = mid, y = after_scale(.85), label = yr),
#       data = dt[yshow == TRUE], angle = 90,
#       inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
#       size = font.size, hjust = "right"
#     )
#   )
# }