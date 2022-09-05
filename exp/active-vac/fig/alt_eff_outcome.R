
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_eff_deaths.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

tar <- match.arg(
  tail(.args, 1) |> basename() |> gsub(pattern = "^.+_(.+)\\..+$", replacement = "\\1"),
  c("inf", "sev", "deaths")
)

#' comes key'd
eff.dt <- readRDS(.args[2])[
  date > "2020-12-01"
][
  outcome == tar
][, .(
  scenario, realization, date, outcome, c.effectiveness
)]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

rm(eff.dt)
gc()

plt.dt[, talloc := factor(
  fifelse(
    pas_alloc == "none", as.character(act_alloc), fifelse(
    pas_alloc == "FL", "FL+ring", as.character(pas_alloc)
  ))
)][, qfac := factor(c("No Q", "Q")[quar+1]) ]

p <- ggplot(plt.dt) + aes(
  x = date, y = c.effectiveness,
  color = act_vac,
  linetype = factor(c("unconditional", "conditional")[inf_con+1]),
  sample = realization
) +
  facet_nested(rows = vars(qfac), cols = vars(talloc)) +
  geom_month_background(
    plt.dt, by = c("qfac", "talloc"),
    font.size = 3, value.col = "c.effectiveness",
    ymax = plt.dt[, max(c.effectiveness)],
    ymin = plt.dt[, min(c.effectiveness)]
  ) +
  stat_spaghetti(
    aes(alpha = after_stat(sampleN^-1))
  ) +
  geom_hline(aes(yintercept=0, color = "none")) +
  scale_y_continuous(
    name = sprintf(
      "Cumulative Effectiveness\nAgainst Incidence of %s",
      switch(tar, inf = "Infection", sev = "Severe Disease", deaths = "Death", stop())
    )
  ) +
  scale_x_null() +
  scale_color_discrete("Active Vax.") +
  scale_linetype_discrete("Conditional Vax.") +
  scale_alpha(range = c(0.01, 1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
