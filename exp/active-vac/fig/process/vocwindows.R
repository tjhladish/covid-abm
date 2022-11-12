
.pkgs <- c("data.table", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args = commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", "vocpattern.rds"),
  file.path("fig", "process", "vocwindows.rds")
))

load(.args[1])

voc.dt <- readRDS(.args[2])

fwd.dt <- voc.dt[
  measure != "vocprevWT", .(
    start = c(
      date[which.max(0.5 < value)],
      date[which.max(0.75 < value)],
      date[which.max(0.95 < value)]
    ),
    q = c(0.5, 0.75, 0.95)
  ),
  by=.(measure, realization)
]

rev.dt <- voc.dt[measure != "vocprevWT"][
  order(date, decreasing = TRUE), .(
    end = c(
      date[which.max(0.5 < value)],
      date[which.max(0.75 < value)],
      date[which.max(0.95 < value)]
    ),
    q = c(0.5, 0.75, 0.95)
  ),
  by=.(measure, realization)
]

fwd.dt[
  rev.dt, on=.(realization, measure, q)
][,.(start = median(start)+1, end = median(end)-1), by=.(measure, q)][,
  mids := start + (end - start)/2
] |>
  store(args = .args, obj = _)
