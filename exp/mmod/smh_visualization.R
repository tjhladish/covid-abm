
# working / invocation directory should be exp/mmod

.args <- if (interactive()) c(
  "smhdata.rds",
  "ccdata.rds",
  "curdata.rds",
  "smccomparison.png"
) else commandArgs(trailingOnly = TRUE)

library(data.table)
library(ggplot2)

# the submitted SMH outputs
smh_dt <- readRDS(.args[1])
smh_dt$type <- NULL
# the gold standard data
cc_dt <- readRDS(.args[2])[, scenario := fifelse(target == "exdeath", "excessdeaths", "reported")][, quantile := 0.5]
cc_dt[target == "exdeath", target := "death"]
# our updated model validation runs
update_dt <- readRDS(.args[3])[,.(
  quantile = c(0.025, 0.5, .975),
  value = quantile(value, probs = c(0.025, 0.5, 0.975))
), by=.(target, date) ][, scenario := "model-rev"]

cc_dt <- cc_dt[between(date, min(update_dt[, min(date)], smh_dt[, min(date)]), max(update_dt[, max(date)], smh_dt[, max(date)])+1)]

p <- ggplot(rbindlist(list(
  smh_dt, cc_dt, update_dt
), use.names = TRUE, fill = TRUE)) + aes(x = date, color = scenario) +
  facet_grid(. ~ target, labeller = labeller(
    target = c(case = "Cases", hosp = "Hospitalizations", death = "Deaths")
  )) +
  geom_ribbon(
    aes(ymin = `0.025`, ymax = `0.975`, color = NULL, fill = scenario),
    data = \(dt) dt[!is.na(quantile) & (quantile != 0.5)] |> dcast.data.table(
      scenario + target + date ~ quantile, value.var = "value"
    ),
    alpha = 0.2
  ) +
#  geom_line(aes(y = value, group = realization), data = \(dt) dt[is.na(quantile) & scenario == "model-rev"], alpha = 1/100) +
  geom_line(aes(y = value), data = \(dt) dt[quantile == 0.5 & !(scenario %in% c("reported", "excessdeaths"))]) +
  geom_point(aes(y = value), data = \(dt) dt[scenario == "reported"], shape = 21) +
  geom_point(aes(y = value), data = \(dt) dt[scenario == "excessdeaths"], shape = 19) +
  theme_minimal() +
  theme(
    legend.position = c(1, 1), legend.justification = c(1, 1),
    legend.text = element_text(margin = margin(t=5, b=5)),
    axis.text.x = element_text(hjust = 0)
  ) +
  scale_color_manual(
    NULL, label = c(
      reported = "Reported",
      excessdeaths = "Excess Deaths",
      "model-rev" = "2023 Model: low CFR, moderate immunity,\nhigh transmissibilty",
      "optSev_highIE" = "UF SMH: low CFR,\nlow immunity & transmissibility",
      "optSev_lowIE" = "UF SMH: low CFR,\nhigh immunity & transmissibility",
      "pessSev_highIE" = "UF SMH: high CFR,\nlow immunity & transmissibility",
      "pessSev_lowIE" = "UF SMH: high CFR,\nhigh immunity & transmissibility"
    ),
    breaks = c("reported", "excessdeaths", "model-rev", "pessSev_lowIE", "optSev_lowIE", "pessSev_highIE", "optSev_highIE"),
    values = c(
      reported = "black", excessdeaths = "black",
      "model-rev" = "dodgerblue",
      "pessSev_lowIE" = "firebrick",
      "optSev_lowIE" = "forestgreen",
      "pessSev_highIE" = "salmon",
      "optSev_highIE" = "yellowgreen"
    ),
    aesthetics = c("color", "fill"),
    guide = guide_legend(override.aes = list(
      shape = c(21, 19, rep(NA, 5)),
      lty = c("blank", "blank", rep("solid", 5)),
      fill = c(NA, NA, "dodgerblue", "firebrick", "forestgreen", "red", "green")
    ))
  ) + scale_x_date(name = NULL) +
  scale_y_continuous(
    "Weekly Incidence per 10K",
    trans = "log10", labels = \(br) sprintf("10^%i", log10(br))
  )

ggsave(tail(.args, 1), p, width = 10, height = 6, bg = "white")