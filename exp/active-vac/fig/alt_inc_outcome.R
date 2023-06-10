
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_inc_inf.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

tar <- match.arg(
  tail(.args, 1) |> basename() |> gsub(pattern = "^.+_(.+)\\..+$", replacement = "\\1"),
  c("inf", "sev", "deaths", "doses")
)

#' comes key'd
eff.dt <- readRDS(.args[2])[
  eval(datefilter)
][
  outcome == tar
][, .(
  scenario, realization, date, outcome, value, averted
)]

eff.dt[order(date),
  c("c.value", "i.c.value") := .(cumsum(value), cumsum(value + averted)),
  by=.(scenario, realization, outcome)
]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

dt2 <- plt.dt[, .(
  act_vac="none", c.value = i.c.value[1]
), by=setdiff(key(plt.dt), "act_vac") ]
setkeyv(dt2, key(plt.dt))

rm(eff.dt)
gc()

plt.qs <- rbind(quantile(
  plt.dt,
  j = .(c.value), sampleby = "realization",
  probs = qprobs(c(`90`=.9), mid = TRUE, extent = FALSE)
  ),
  quantile(
    dt2,
    j = .(c.value), sampleby = "realization",
    probs = qprobs(c(`90`=.9), mid = TRUE, extent = FALSE)
  )
)[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LS", "IS", "HS", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

p <- ggplot(plt.qs) + aes(
  x = date,
  color = act_vac,
  linetype = factor(c("unconditional", "conditional")[inf_con+1])
) +
  facet_nested(rows = vars(qfac), cols = vars(talloc)) +
  geom_month_background(
    plt.qs, by = c("qfac", "talloc"),
    font.size = 3, value.col = "qmed",
    ymax = plt.qs[, max(q90h)],
    ymin = plt.qs[, min(q90l)]
  ) +
  geom_ribbon(aes(ymin=q90l, ymax=q90h, fill=act_vac, color=NULL), alpha=0.25) +
  geom_line(aes(y=qmed)) +
  scale_color_discrete(
    "Vaccine Program",
    breaks = c("none", "risk", "ring"),
    labels = c(
      none = "Standard Program",
      ring="Ring Vaccination",
      risk="Risk-Based Strategy"
    ),
    aesthetics = c("color", "fill")
  ) +
  scale_y_continuous(
    name = sprintf(
      "Per 10k, Cumulative\nIncidence of %s",
      switch(tar, inf = "Infection", sev = "Severe Disease", deaths = "Death", doses = "Vaccination", stop())
    )
  ) +
  scale_x_null() +
  scale_linetype_manual("Vaccinate...", labels = c(conditional="Condtional on\nCase History", unconditional = "Unconditionally"), values = c(conditional="dashed", unconditional="solid")) +
  scale_alpha(range = c(0.02, 1)) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1), strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
