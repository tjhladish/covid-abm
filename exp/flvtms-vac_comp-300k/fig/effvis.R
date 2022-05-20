
.pkgs <- c("data.table", "ggplot2")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args <- if (interactive()) c(
  file.path("tmp", "digest.rds"),
  file.path("fig", "effectiveness.png")
) else commandArgs(trailingOnly = TRUE)

dt <- readRDS(.args[1])

non.dt <- dt[vac == 0, .SD, .SDcols = c("realization", "variable", "week", "value")]
int.dt <- dt[vac != 0, .SD, .SDcols = c("realization", "state", "variable", "week", "value")]

eff.dt <- int.dt[non.dt, on=.(realization, variable, week)]

eff.dt[order(week), cv := cumsum(i.value) - cumsum(value), by=.(realization, state, variable)]

#' TODO decide appropriate week to measure from
intweek <- 43

eff.dt[week > intweek, eff := 1 - cumsum(value)/cumsum(i.value), by=.(realization, state, variable) ]

varkey <- c(c="Symptomatic Infections", s="Severe Infections", d="Deaths")

eff.dt[, varf := factor(varkey[variable], levels = varkey, ordered = TRUE)]

obs.p <- ggplot(eff.dt) + aes(week, value, color = state, group = interaction(state, realization)) +
  facet_grid(varf ~ ., scale = "free_y") +
  geom_line(alpha = 0.1) +
  geom_line(mapping = aes(y=i.value, color = "NONE"), data = function(dt) dt[state == "FL"], alpha = 0.1) +
  theme_minimal() +
  scale_y_continuous("Incidence") +
  scale_color_discrete(name = NULL, guide = guide_legend(override.aes = list(alpha = 1))) +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1)
  )

averted.p <- ggplot(eff.dt) + aes(week, cv/1000, color = state, group = interaction(state, realization)) +
  facet_grid(varf ~ ., scale = "free_y") +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  scale_y_continuous("Cumulative 1K averted") +
  scale_color_discrete(name = NULL, guide = guide_legend(override.aes = list(alpha = 1))) +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1)
  )

eff.p <- ggplot(eff.dt[week > intweek]) + aes(week, eff, color = state, group = interaction(state, realization)) +
  facet_grid(varf ~ .) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  scale_y_continuous("Cumulative Fraction Averted") +
  scale_color_discrete(name = NULL, guide = guide_legend(override.aes = list(alpha = 1))) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1)
  )

#' TODO refactor these to independent scripts
ggsave("incidence.png", obs.p, width = 6, height = 8, units = "in", dpi = 600, bg = "white")
ggsave("cum_averted.png", averted.p, width = 6, height = 8, units = "in", dpi = 600, bg = "white")
ggsave("cum_effectiveness.png", eff.p, width = 6, height = 8, units = "in", dpi = 600, bg = "white")
