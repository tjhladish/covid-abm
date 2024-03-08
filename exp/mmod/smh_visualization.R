
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
# our updated model outputs
update_dt <- readRDS(.args[3])[, scenario := "model-rev"][, quantile := 0.5]

cc_dt <- cc_dt[date <= max(update_dt[, max(date)], smh_dt[, max(date)])]

p <- ggplot(rbindlist(list(
  smh_dt, cc_dt, update_dt
), use.names = TRUE)) + aes(x = date, color = scenario) +
  facet_grid(. ~ target, labeller = labeller(
    target = c(case = "Cases", hosp = "Hospitalizations", death = "Deaths")
  )) +
  geom_ribbon(
    aes(ymin = `0.025`, ymax = `0.975`, color = NULL, fill = scenario),
    data = \(dt) dt[quantile != 0.5] |> dcast.data.table(
      scenario + target + date ~ quantile, value.var = "value"
    ),
    alpha = 0.2
  ) +
  geom_line(aes(y = value), data = \(dt) dt[quantile == 0.5 & !(scenario %in% c("reported", "excessdeaths"))]) +
  geom_point(aes(y = value), data = \(dt) dt[scenario %in% c("reported", "excessdeaths")]) +
  theme_minimal() +
  theme(
    legend.position = c(1, 1), legend.justification = c(1, 1)
  ) +
  scale_color_discrete(
    "Scenario", label = c(
      reported = "Reported",
      excessdeaths = "Excess Deaths",
      "model-rev" = "Ultimate Model",
      "optSev_highIE" = "SMH, low Severity, high Imm. Escape",
      "optSev_lowIE" = "SMH, low Severity, low Imm. Escape",
      "pessSev_highIE" = "SMH, high Severity, high Imm. Escape",
      "pessSev_lowIE" = "SMH, high Severity, low Imm. Escape"
    ), aesthetics = c("color", "fill")) +
  scale_x_date(name = NULL) +
  scale_y_continuous(
    "Incidence per 10K",
    trans = "log10", labels = \(br) sprintf("10^%i", log10(br))
  )

ggsave(tail(.args, 1), p, width = 8, height = 6, bg = "white")