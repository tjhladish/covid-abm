
library(data.table)
library(ggplot2)

.args <- if (interactive()) c(
  "round11_detection_2020-02-05_start.txt",
  "current_detection_2020-02-10_start.txt",
  "smh_detection_comparison.png"
) else commandArgs(trailingOnly = TRUE)

read_det_file <- function(pth) {
  date_ref <- as.Date(gsub(".*_(2020.*)_start.*", "\\1", pth))
  dt <- fread(pth)[, date := date_ref + V1]
}

smh_dt <- read_det_file(.args[1])
cur_dt <- read_det_file(.args[2])

