
.pkgs <- c("data.table", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args = commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", "vocpattern.rds"),
  file.path("fig", "process", "vocwindows.rds")
))

load(.args[1])

voc.dt <- readRDS(.args[2])

voc.dt[measure != "vocprevWT"][, .(
  s50 = date[which.max(value > 0.25)],
  e50 = date[which.max(value > 0.75)],
  s90 = date[which.max(value > 0.05)],
  e90 = date[which.max(value > 0.95)],
  s95 = date[which.max(value > 0.025)],
  e95 = date[which.max(value > 0.975)]
), by=.(measure, realization)][, .(
  start = c(
    mean(s50), mean(s90), mean(s95)
  ),
  end = c(
    mean(e50), mean(e90), mean(e95)
  ),
  mids = c(
    mean(s50 + (e50-s50)/2), mean(s90 + (e90-s90)/2), mean(s95+ (e95-s95)/2)
  ),
  q = c(.5, .9, .95)
), by = .(measure)
] |> store(args = .args, obj = _)
