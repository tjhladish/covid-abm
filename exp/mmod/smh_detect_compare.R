
library(data.table)
library(ggplot2)
library(geomtextpath)
library(scales)
library(patchwork)

.args <- if (interactive()) c(
  "round11_detection_2020-02-05_start.txt",
  "current_detection_2020-02-10_start.txt",
  "smh_reporting.png"
) else commandArgs(trailingOnly = TRUE)

read_det_file <- function(pth) {
  date_ref <- as.Date(gsub(".*_(2020.*)_start.*", "\\1", pth))
  dt <- fread(
    pth, col.names = c("day", "asymp", "mild", "severe", "crit", "death")
  )[, date := date_ref + day][,
       mild := asymp + mild - asymp*mild
  ][,
       severe := mild + severe - mild*severe
  ][,
       crit := severe + crit - severe*crit
  ][,
       death := crit + death - crit*death
  ]
  dt[, .(date, asymp, mild, severe, crit, death)]
}

prepare <- function(..., id.vars = c("realization", "date")) setkeyv(melt(
  rbind(..., fill = TRUE),
  id.vars = id.vars,
  variable.name = "measure", variable.factor = FALSE
), c("measure", id.vars))

smh_dt <- prepare(read_det_file(.args[1])[, .(
  realization = 1, date, asymp, mild, severe, crit
)])
cur_dt <- prepare(read_det_file(.args[2])[, .(
  realization = 1, date, asymp, mild, severe, crit
)])

setnames(smh_dt, "measure", "outcome")
smh_dt[, measure := "detection"]
setnames(cur_dt, "measure", "outcome")
cur_dt[, measure := "detection"]

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

geom_month_bars <- function(
    data,
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    ytrans = "identity",
    datafn = tsref.dt,
    ...
) {
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]
  
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
    )
  )
}

geom_grid <- function(
    ys, FUN = as.character
) {
  geom_texthline(
    aes(yintercept = ys, label = label),
    data = data.table(ys = ys, label = FUN(ys)), inherit.aes = FALSE,
    gap = TRUE, hjust = 0.325, color = "grey65", fontface = "bold"
  )
}

p.core <- function(
    dt, ymin = NA, ymax = NA, ytrans = "identity", gridy
) ggplot(dt) +
  aes(date, value) +
  geom_grid(gridy) +
  geom_month_bars(dt, ymin = ymin, ymax = ymax, ytrans = ytrans) +
  theme_minimal() + theme(text = element_text(face = "bold"), axis.title.y = element_text(size = rel(1.5))) +
  scale_x_null()

match.call.defaults <- function(
    definition = sys.function(sys.parent()),
    call = sys.call(sys.parent()),
    expand.dots = TRUE,
    envir = parent.frame(2L)
) {
  # get the call
  mc <- match.call(definition, call, expand.dots, envir)
  # get the formals, tossing any ellipsis
  fs <- formals(definition, envir)
  fs$... <- NULL
  
  # for any arguments set in formals & not in the call
  for(nm in setdiff(names(fs), names(mc)))
    mc[nm] <- fs[nm] # add those to the call
  
  return(mc)
}

rejig <- function(FUN, ..., .ENV = environment(FUN)) {
  # initial validation
  stopifnot(
    "FUN isn't a function." = is.function(FUN),
    "FUN is a primitive function." = !is.primitive(FUN)
  )
  
  dots <- as.list(match.call())[-1] # get some new defaults
  dots$FUN <- dots$.ENV <- NULL # drop all the not-defaults
  
  if (length(dots) == 0) {
    warning("... is empty. Just returning FUN.")
    return(FUN)
  }
  
  .FUN <- FUN # make a duplicate of FUN
  forms <- formals(FUN) # get the original defaults
  
  # potentially more validation: check for ... argument
  # in FUN and try to partial match all arguments in
  # rejig
  hasdots <- "..." %in% names(forms)
  replacements <- names(forms)[pmatch(names(dots), names(forms))]
  
  if (any(is.na(replacements)) && !hasdots) {
    errmsg <- sprintf("
FUN does not have ... argument, and
rejig ... arguments do not match FUN arguments:
%s
", names(dots)[is.na(replacements)] |> paste(collapse = ", ")
    )
    stop(errmsg)
  }
  
  # correct any partially matched defaults
  names(dots)[!is.na(replacements)] <- replacements[!is.na(replacements)]
  # set the new defaults
  formals(.FUN)[names(dots)] <- dots
  environment(.FUN) <- .ENV
  
  if (hasdots && any(is.na(replacements))) {
    # the internals of FUN may pass around the ellipsis, which now
    # excludes newly set default variables, so need to use it
    body(.FUN) <- substitute({
      mc <- cabputils::match.call.defaults()
      mc[[1]] <- FUN
      eval(mc)
    })
  }
  
  return(.FUN)
  
}

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

p.detect_smh <- p.core(
  smh_dt[date < "2022-04-01"], ymin = 0, ymax = 1, gridy = c(0.25, 0.5, 0.75)
) + aes(linetype = outcome) +
  geom_line() +
  scale_y_fraction(name = "UF SMH Round 11\nProbability of Detection", guide = "none") +
  geom_text(mapping = aes(label=lab, linetype = NULL, shape = NULL, hjust = align), data = data.table(
    lab = c("Asymptomatic", "Mild", "Severe", "Critical"),
    date = as.Date(c("2020-07-15", "2020-05-10", "2020-04-15", "2020-08-22")),
    value = c(.1, .22, .65, .9),
    align = c(0.5, 1, 0.5, 0.5),
    measure = "detection"
  ), size = 5) +
  scale_linetype_manual(
    name = NULL, breaks = rev(c("asymp", "mild", "severe", "crit")),
    labels = c(asymp = "Asymptomic", mild = "Mild", severe = "Severe", crit = "Critical"),
    values = c(asymp = "dotted", mild = "dotdash", severe = "longdash", crit = "solid")
  ) +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank()
  ) + coord_cartesian(clip = "off")

p.detect_cur <- p.core(
  cur_dt[date < "2022-04-01"], ymin = 0, ymax = 1, gridy = c(0.25, 0.5, 0.75)
) + aes(linetype = outcome) +
  geom_line() +
  scale_y_fraction(name = "Ultimate Model\nProbability of Detection", guide = "none") +
  geom_text(mapping = aes(label=lab, linetype = NULL, shape = NULL, hjust = align), data = data.table(
    lab = c("Asymptomatic", "Mild", "Severe", "Critical"),
    date = as.Date(c("2020-07-15", "2020-05-10", "2020-04-15", "2020-08-22")),
    value = c(.18, .2, .675, .9),
    align = c(0.5, 1, 0.5, 0.5),
    measure = "detection"
  ), size = 5) +
  scale_linetype_manual(
    name = NULL, breaks = rev(c("asymp", "mild", "severe", "crit")),
    labels = c(asymp = "Asymptomic", mild = "Mild", severe = "Severe", crit = "Critical"),
    values = c(asymp = "dotted", mild = "dotdash", severe = "longdash", crit = "solid")
  ) +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank()
  ) + coord_cartesian(clip = "off")

geom_month_end <- function(
    data,
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    m.labels = m.abb,
    font.size = 8, font.face = "bold",
    ytrans = "identity",
    datafn = tsref.dt,
    m.y = 0.0, y.y = 0.3,
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
      size = font.size,
      vjust = if (y.y > 0) "bottom" else "top",
      fontface = font.face
    ),
    geom_text(
      mapping = aes(x = start+5, y = xform(ymax, ymin, y.y), label = yr),
      data = dt[yshow == TRUE],
      inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
      size = font.size, hjust = 0, fontface = font.face,
      vjust = if (y.y > 0) "top" else "bottom"
    )
  )
}

conserved <- list(
  scale_x_null(),
  theme_minimal(),
  theme(text = element_text(face = "bold"))
)

m.abb <- gsub("^(.).+$","\\1", month.abb)

p.bot <- ggplot() +
  geom_month_end(
    smh_dt[date < "2022-04-01"][, .(date, value)],
    ymax = 1, ymin = 0, m.y = -0.05, y.y = -0.6
  ) + conserved + theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank()
  ) + coord_cartesian(ylim=c(-0.6,0))

p.fin <- (p.detect_smh / p.bot / p.detect_cur) + plot_layout(nrow = 3, heights = c(1, .2, 1))

ggsave(tail(.args, 1), p.fin, height = 5*2, width = 8*2, bg = "white")
