
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils", "geomtextpath")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("alt_eff.rds", "digest-key.rds")),
  file.path("fig", "output", "alt_eff_all.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  eval(datefilter)
][
  outcome %in% c("inf", "sev", "deaths")
][, .(
  scenario, realization, date, outcome, c.effectiveness
)]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])[scenario %in% intscns]

plt.dt <- setkeyv(
  eff.dt[scn.dt, on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)[inf_con == FALSE]

rm(eff.dt)
gc()

plt.dt[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LIC", "MIC", "HIC", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

plt.qs <- quantile(
  plt.dt,
  j = .(c.effectiveness), sampleby = "realization",
  probs = qprobs(c(`90`=.9), mid = TRUE, extent = FALSE)
)[, talloc := factor(
  fifelse(
    pas_alloc == "none",
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LIC", "MIC", "HIC", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

p <- ggplot(plt.qs) + aes(
  x = date,
  color = act_vac,
  linetype = factor(c("nonpi", "wquar")[quar+1])
) +
  facet_nested(
    rows = vars(outcome), cols = vars(talloc), switch = "y",
    scales = "free_y", labeller = labeller(
      outcome = c(inf = "Infection", sev="Severe Disease", deaths = "Deaths")
    )
  ) +
  geom_month_background(
    plt.qs, by = c(row="outcome", col="talloc"),
    font.size = 3, value.col = "qmed", max.col = "q90h", min.col = "q90l"
  ) +
  geom_ribbon(aes(ymin=q90l, ymax=q90h, fill=act_vac, color=NULL), alpha=0.25) +
  geom_line(aes(y=qmed)) +
  scale_color_discrete(
    "Vaccine Program",
    breaks = c("risk", "ring"),
    labels = c(
      ring="Ring Vaccination",
      risk="Risk-Based Strategy"
    ),
    aesthetics = c("color", "fill")
  ) +
  geom_hline(
    aes(yintercept = 0, color = "none"),
    data = \(dt) dt[,.SD[1],by=.(outcome, talloc)]
  ) +
  geom_texthline(
    aes(yintercept = 0, color = "none", label = "Standard\nProgram"),
    inherit.aes = FALSE, show.legend = FALSE, data = \(dt) dt[outcome == "sev"][talloc == "LIC",.SD[1],by=.(outcome, talloc)],
    hjust = 0, gap = FALSE
  ) +
  scale_y_continuous(
    name = "Cumulative Effectiveness\nAgainst Incidence of ..."
  ) +
  scale_x_null() +
  scale_linetype_manual(
    "Extra NPI", labels = c(nonpi="None", wquar = "Quarantine Contacts"),
    values = c(nonpi="dashed", wquar="solid")
  ) +
  scale_alpha(range = c(0.01, 1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
