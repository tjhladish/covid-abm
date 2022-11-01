
.pkgs <- c("data.table", "ggplot2", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args <- commandArgs(args = c("~/Downloads", "something.png"))

dt <- .args[1] |> list.files(pattern = "vpdr", full.names = TRUE) |>
  lapply(fread) |> rbindlist() |> melt.data.table(id.vars = c("serial", "day"))

dt[, c("dose", "var") := tstrsplit(variable, split = "_") ]

p <- ggplot(dt[day > 350]) + aes(x=day + as.Date("2020-02-10"), y=value, color=dose) +
  facet_grid(var ~ serial, scales = "free_y", switch = "y") +
  geom_line() +
  theme_minimal() + theme(
    legend.position = "bottom", strip.text = "outside"
  ) +
  scale_x_date(
    name = NULL, date_breaks = "3 months", date_labels = "%b %y", minor_breaks = NULL
  )

store(.args, p, width = 10, height = 6, bg = "white")
