
.pkgs <- c("data.table", "ggplot2")
stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args <- if (interactive()) c(
  "covax_doses.txt",
  "covax_dose_compare.png"
) else commandArgs(trailingOnly = TRUE)

dt <- fread(.args[1]) |>
  melt(id.vars = "date") |>
  (\(dt) dt[order(date), cvalue := cumsum(value), by=variable])() |>
  melt(id.vars = c("date","variable"), variable.name = "measure")

p <- ggplot(dt) + aes(
  x = date, y = value, color = variable
) + facet_grid(
  measure ~ ., scale = "free_y", switch = "y",
  labeller = labeller(measure = c(value = "Incident", cvalue = "Cumulative"))
) +
  geom_step(data = \(dt) dt |> subset(variable != "provisioned")) +
  geom_point(data = \(dt) dt |> subset(variable == "provisioned"), alpha = 0.2) +
  theme_minimal() +
  theme(
    legend.position = c(0.1, .9), legend.justification = c(0, 1),
    strip.placement = "outside"
  ) +
  scale_x_date(name = NULL, date_breaks = "months", date_labels = "%b %C", minor_breaks = NULL) +
  scale_y_continuous("Per 10k, age 5+") +
  scale_color_manual(
    name = NULL, breaks = c(
      "assumeddoses", "assumedcourses", "provisioned"
    ),
    labels = c(
      assumedcourses = "Sindh Analysis (Courses)",
      assumeddoses = "Sindh Analysis (Doses)",
      provisioned = "UNICEF"
    ),
    guide = guide_legend(override.aes = list(
      linetype=c("solid", "solid", "blank"), shape = c(NA, NA, 20), alpha = 1
    )),
    values = c(assumeddoses = "dodgerblue", assumedcourses = "firebrick", provisioned = "black")
  )

ggsave(tail(.args, 1), p, width = 10, height = 5, bg = "white")
