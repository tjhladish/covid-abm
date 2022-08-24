.pkgs <- c("data.table")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args = if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  "rcasedeath-florida.csv",
  file.path("fig", "input", "outcomes.rds")
) else commandArgs(trailingOnly=TRUE)

load(.args[1])

tarpop <- gsub("^.+-(\\w+)\\.csv$", "\\1", .args[2])
per10k <- 1e4/pop[tarpop]

ed <- fread(.args[2])[Date <= endday] |> setnames("Date", "date")
ed[, rcase := rcase * per10k ]
#' TODO fix upstream to only report on wday == 7? if so, need to change 7 below
ed[, rdeath := excess*per10k*7 ]
ed[wday(date) != 7, rdeath := NA ]
ed[, crcase := cumsum(rcase) ]
ed[!is.na(rdeath), crdeath := cumsum(rdeath) ]

saveRDS(ed, tail(.args, 1))
