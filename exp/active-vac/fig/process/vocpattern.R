
.pkgs <- c("data.table", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args = commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", "validation.rds"),
  file.path("fig", "process", "vocpattern.rds")
))

load(.args[1])

.args[2] |> readRDS() |> DT(,.(
  realization, date, vocprev1, vocprev2, vocprev3,
  vocprevWT = pmax(1-(vocprev1 + vocprev2 + vocprev3), 0)
)) |> prepare() |> store(.args, obj = _)
