
.args <- if (interactive()) c("plottools.rda") else commandArgs(trailingOnly = TRUE)

#' `gg.scale.wrapper` enables easy creation of re-usable ggplot scales
#' (a Facade Factory, if you're a fan of Go4 Design Patterns)
#' 
#' @param scale_fun the ggplot scale function that will ultimately be called
#' @param ... whatever default arguments for that scale function; the
#'   ones you define will override the defaults
gg.scale.wrapper <- function(
  scale_fun,
  ...
) {
  stopifnot(!missing(scale_fun))
  defs <- list(...)
  if (!length(defs)) warning(
    "provided no default arguments; consider using scale_fun directly."
  )
  
  return(function(...) {
    #' this different ... is how you get a function back that let's you
    #' override defaults, set other arguments to scale_... functions
    .ellipsis <- list(...)
    .args <- defs
    .args[names(.ellipsis)] <- .ellipsis
    do.call(scale_fun, .args)
  })
}

#' then use like this to create your scale-with-defaults:
scale_x_Ms <- gg.scale.wrapper(
  scale_x_date,
  name = NULL,
  breaks = function(l) {
    rl <- round(as.POSIXlt(l), "month")
    as.Date(seq(rl[1] + 60*60*24*15, rl[2], by="month"))
  },
  labels = function(bs) {
    gsub("^(.).+$","\\1",month.abb[month(bs)])
  },
  minor_breaks = NULL
)

month_background_scale <- function(
  ylims = c(0,Inf),
  band.color = "grey",
  band.alpha = 0.25,
  ...
) list(
  geom_rect(
    aes(xmin = start, xmax=end, ymax=ylims[2], ymin=ylims[1]),
    data = function(dt) {
      dr <- round(as.POSIXlt(dt[, range(date)]), "month")
      pts <- seq(dr[1], dr[2], by="month")
      if (length(pts) %% 2) pts <- head(pts, -1)
      data.table(
        start = as.Date(pts[seq(1, length(pts), by=2)]),
        end   = as.Date(pts[seq(2, length(pts), by=2)])
      )
    },
    fill = band.color, alpha = band.alpha,
    inherit.aes = FALSE, show.legend = FALSE
  ),
  scale_x_Ms(...)
)

event_vline <- function(
  vline.date = as.Date("2020-01-01"),
  vline.lty = "dotted",
  vline.color = "dodgerblue",
  labelf = function(d) strftime(d, "Vaccination Start,\n%Y-%m-%d"),
  dataf = function(dt) dt[
    variable == "R[t]", .(
      intervention = TRUE,
      date = vline.date,
      value = max(value)
    ), by = .(variable, var)
  ]
) list(
  geom_vline(
    xintercept = vline.date, linetype = vline.lty, color = vline.color
  ),
  if (!is.null(labelf)) geom_text(
    aes(label = labelf(date)),
    data = dataf,
    show.legend = FALSE,
    hjust = "left",
    vjust = "top",
    size = 2, alpha = 0.5,
    nudge_x = 7
  )
)

scale_color_intervention <- gg.scale.wrapper(
  scale_color_manual,
  name = NULL,
  values = c(`FALSE`="black", `TRUE`="dodgerblue"),
  labels = c(`FALSE`="No Vaccine", `TRUE`="Vaccine")
)

save(list=ls(), file = tail(.args, 1))