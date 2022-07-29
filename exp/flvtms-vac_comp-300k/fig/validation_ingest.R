
.pkgs <- c("data.table")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "plotlog"),
  file.path("fig", "validation.rds")
) else commandArgs(trailingOnly = TRUE)

fls <- list.files(.args[1], full.names = TRUE)

extractor <- function(fl) {
  ret <- setnames(fread(fl)[,
    c("crcase", "crdeath") := .(cumsum(rcase), cumsum(rdeath))
  ][,
    VESAvg := filter(VES, rep(1/7, 7), sides = 2)
  ][,
    brkthruRatioAvg := filter(brkthruRatio, rep(1/7, 7), sides = 2)
  ][,
    c("tot_std_doses", "tot_urg_doses") := .(cumsum(std_doses), cumsum(urg_doses))
  ], "serial", "realization")

  ret[wday(date) != 7, crdeath := NA]
  ret[, rdeath := NA ]
  ret[wday(date) == 7, rdeath := c(crdeath[1], diff(crdeath)) ]
  ret
}

dt <- rbindlist(lapply(fls, extractor))

saveRDS(dt, tail(.args, 1))

