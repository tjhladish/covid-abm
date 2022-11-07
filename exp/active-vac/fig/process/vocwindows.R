
.pkgs <- c("data.table", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args = commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", "vocpattern.rds"),
  file.path("fig", "process", "vocwindows.rds")
))

load(.args[1])

.args[2] |> readRDS() |> store(args, obj = _)
