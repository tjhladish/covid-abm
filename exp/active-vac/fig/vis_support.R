
#' DANGER
#' because this ends with saving everything, have to clean env first when
#' used interactively.
#' TODO: do via creating and saving an environment, rather than everything
rm(list = ls())

.pkgs <- c("data.table", "ggplot2", "ggrepel", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args <- if (interactive()) {
  file.path("fig", "vis_support.rda")
} else commandArgs(trailingOnly = TRUE)

# MAGIC DATE
endday <- as.Date("2022-03-07")
vendday <- as.Date("2022-03-31")
startday <- as.Date("2020-12-01")
datefilter <- expression(between(date, startday, endday))
outfilter <- expression(outcome %in% c("inf", "deaths"))
pop <- c(florida = 21538187, escambia = 312212, dade = 2794464)

scale_y_effectiveness <- rejig(
  scale_y_continuous,
  name = "Cumulative Effectiveness Against ...",
  breaks = seq(0, 1, by=0.25)
)

scale_y_doses <- rejig(
  scale_y_continuous,
  name = "Per 10k,\nDaily Doses Administered"
)

scale_y_cdoses <- rejig(
  scale_y_continuous,
  name = "Per 10k,\nCumulative Doses Administered"
)

scale_y_incidence <- rejig(
  scale_y_continuous,
  name = "Per 10k,\nIncidence of ..."
)

scale_y_averted <- rejig(
  scale_y_continuous,
  name = "Per 10k,\nAverted Incidence of ..."
)

scale_x_null <- rejig(
  scale_x_continuous,
  name = NULL, breaks = NULL, minor_breaks = NULL,
  expand = expansion(
    mult = c(0, 0), add = c(0, 0)
  )
)

scale_y_fraction <- rejig(
  scale_y_continuous,
  name = "Fraction of Population",
  breaks = (0:4)/4,
  limits = c(0, 1),
  expand = expansion(
    mult = c(0, 0), add = c(0, 0)
  )
)

scale_linetype_quar <- rejig(
  scale_linetype_manual,
  name = "Extra NPI", labels = c(nonpi="None", wquar = "Quarantine Contacts"),
  values = c(wquar="dashed", nonpi="solid"),
  guide = guide_legend(title.position = "top", title.hjust = 0.5)
)

scale_color_strategy <- rejig(
  scale_color_manual,
  name = "Vaccine Program",
  breaks = c("none", "risk", "age", "ring"),
  labels = c(
    none="Standard Program",
    ring="Ring Vaccination",
    age = "Age-Based Strategy",
    risk="Risk-Based Strategy"
  ),
  aesthetics = c("color", "fill"),
  values = c(none = "black", ring = "#fb6502", risk = "#00529b", age = "#006b35"),
  guide = guide_legend(title.position = "top", title.hjust = 0.5)
)

scale_linetype_scenario <- rejig(
  scale_linetype_manual,
  name = "Intervention",
  breaks = c("none", "qonly", "vonly", "vandq"),
  labels = c(vandq="Vac.+Quar.", vonly = "Vac. Only", qonly = "Quar. Only", none = "None"),
  values = c(vandq="solid", vonly = "dashed", qonly = "dotted", none = "solid"),
  drop = TRUE, limits = force,
  guide = guide_legend(override.aes = list(color = "grey"))
)

scale_color_scenario <- rejig(
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

scale_linetype_fullscenario <- rejig(
  scale_linetype_manual,
  name = "Intervention",
  breaks = fullscnbreaks,
  labels = fullscnlabels,
  values = fullscnlines,
  drop = TRUE, limits = force
  , guide = guide_legend(override.aes = list(fill = NA))
)

scale_color_fullscenario <- rejig(
  scale_color_manual,
  name = "Intervention",
  breaks = fullscnbreaks,
  labels = fullscnlabels,
  values = fullscncolors,
  drop = TRUE, limits = force,
  aesthetics = c("color", "fill")
  , guide = guide_legend(override.aes = list(fill = NA))
)

scale_color_measure <- rejig(
  scale_color_manual,
  name = NULL,
  values = c(
    case = "royalblue3", death = "green4", infection = "orangered",
    hospPrev = "orange", hospInc = "orange", vaxHosp = 'dodgerblue4',
    brkthru = "tan4"
  ),
  guide = "none"
)

scale_color_inputs <- rejig(
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

scale_alpha_measure <- rejig(
  scale_alpha_manual,
  name = NULL,
  breaks = c("observed", "sample", "central", "quantile"),
  labels = measlbls,
  values = c(observed = 0.6, sample = 0.05, quantile = 0.5, central = 1),
  drop = TRUE, limits = force,
  guide = guide_legend(
    override.aes = list(
      linetype = c("blank", "solid", "solid", "blank"),
      size = c(10, 1, 3, 0)
    )
  )
)

scale_shape_measure <- rejig(
  scale_shape_manual,
  name = NULL,
  breaks = c("observed", "sample", "central"),
  labels = measlbls[c("observed", "sample", "central")],
  values = c(observed = 20, sample = NA, central = NA),
  guide = guide_legend(
    override.aes = list(
      linetype = c("blank", "solid", "solid"),
      alpha = c(1, 1/20, 1)
    )
  )
)

geom_observation <- rejig(
  geom_point,
  mapping = aes(shape = "observed"),
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
    ymax = NA, ymin = NA, max.col = value.col, min.col = value.col
) {
  dt <- as.data.table(dt)
  maxer <- function(x) min(max(x, ymax, na.rm = TRUE), ymax, na.rm = TRUE)
  miner <- function(x) max(min(x, ymin, na.rm = TRUE), ymin, na.rm = TRUE)
  if(!is.null(by)) { if ("col" %in% names(by)) {
    #' need to consider ymax across rows if by["col"] exists
    ymref <- dt[, .(mx = maxer(get(max.col)), mn = miner(get(min.col))), by = eval(unname(by["row"]))]
    dt[ymref, c("ym", "yn") := .(mx, mn), on=unname(by["row"])]
    dt <- dt[!is.na(get(by["col"]))]
  } else {
    dt[, c("ym","yn") := .(
      maxer(get(max.col)),
      miner(get(min.col))
    ), by = eval(unname(by)) ]
  } } else {
    dt[, c("ym","yn") := .(
      maxer(get(max.col)),
      miner(get(min.col))
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
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    m.labels = m.abb,
    font.size = 8, font.face = "bold",
    ytrans = "identity",
    datafn = tsref.dt,
    m.y = 0.85, y.y = 0.81,
    text.col = rep("darkgrey", length(col.cycle)),
    ...
) {
  names(text.col) <- c(tail(names(col.cycle), -1), names(col.cycle)[1])
  text.col[is.na(text.col)] <- "white"
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]
  dt[, col := text.col[mtype] ]

  xformer <- scales::as.trans(ytrans)

  xform <- function(ymax, ymin, dropto) {
    xformer$inverse(
      (xformer$transform(ymax)-xformer$transform(ymin))*dropto + xformer$transform(ymin)
    )
  }

  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (xformer$name %like% "log") 0 else -Inf, ymax = if (xformer$name %like% "prob") 1 else Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    ),
    geom_text(
      mapping = aes(x = mid, y = xform(ymax, ymin, m.y), label = m.abb[mon]),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, color = dt$col,
      size = font.size, vjust = "bottom", fontface = font.face
    ),
    geom_text(
      mapping = aes(x = mid, y = xform(ymax, ymin, y.y), label = yr),
      data = dt[yshow == TRUE], angle = 90,
      inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
      size = font.size, hjust = "right", fontface = font.face
    )
  )
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

facet_typical <- rejig(
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
  mapv$shape <- quo(NULL)
  maph <- mapping
  maph$ymax <- NULL; maph$ymin <- NULL
  maph$shape <- quo(NULL)
  res <- list(
    geom_linerange(mapv, data = data, show.legend = FALSE),
    geom_linerange(maph, data = data, show.legend = FALSE)
  )
  showpoint <- !is.null(mapping$shape)
  if (showpoint) {
    mapp <- mapping
    mapp$xmax <- NULL
    mapp$xmin <- NULL
    mapp$ymax <- NULL
    mapp$ymin <- NULL
    # mapp$alpha <- mapp$shape
    res[[length(res)+1]] <- force(geom_observation(mapping = mapp, data = seroprev))
  }
  return(res)
}

prepare <- function(...) setkey(melt(
  rbind(..., fill = TRUE),
  id.vars = c("realization", "date"),
  variable.name = "measure", variable.factor = FALSE
), measure, realization, date)

allplot <- function(
  data.qs, yl, withRef = FALSE,
  col.breaks = if (withRef) c("risk", "age", "ring") else c("none", "risk", "age", "ring"),
  withBands = NULL
) {
  res <- ggplot(data.qs) + aes(
  x = date,
  color = act_vac,
  linetype = factor(c("nonpi", "wquar")[quar+1])
) +
  facet_nested(
    rows = vars(outcome), cols = vars(talloc), switch = "y",
    scales = "free_y", labeller = labeller(
      outcome = c(inf = "Infection", sev = "Severe Disease", deaths = "Deaths", vaccine = "Per 10K,\nDoses Administered")
    )
  ) +
  geom_month_background(
    data.qs, by = c(row="outcome", col="talloc"),
    font.size = 3, value.col = "qmed", max.col = "q90h", min.col = "q90l"
  )

  if (withRef) {
    res <- res +
      geom_hline(
        aes(yintercept = 0, color = "none", linetype = "nonpi"),
        show.legend = FALSE, data = \(dt) dt[,.SD[1],by=.(outcome, talloc)]
      ) +
      geom_texthline(
        aes(yintercept = 0, color = "none", label = "Reference\nProgram"),
        inherit.aes = FALSE, show.legend = FALSE, data = \(dt) dt[talloc == "LIC",.SD[1],by=.(outcome, talloc)],
        hjust = 0, gap = FALSE
      )
  }

  res <- res + geom_ribbon(aes(ymin=q90l, ymax=q90h, fill=act_vac, color=NULL), alpha=0.15) +
#  geom_ribbon(aes(ymin=q50l, ymax=q50h, fill=act_vac, color=NULL), alpha=0.25) +
  geom_line(aes(y=qmed)) +
  scale_color_strategy() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    name = yl, expand = c(0, 0),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  scale_x_null() +
  scale_linetype_quar() +
#  scale_alpha(range = c(0.02, 1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    legend.direction = "horizontal"
  )
  if (!is.null(withBands)) {
# TODO start-end, era

  }
  return(res)
}

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
