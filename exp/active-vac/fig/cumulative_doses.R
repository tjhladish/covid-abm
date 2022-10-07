
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest-doses.rds", "digest-key.rds", "alt_eff.rds")),
  file.path("fig", "output", "cumulative_doses.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#incovax <- fread(.args[5])

tar <- "doses"

scn.dt <- readRDS(.args[3])
scns <- scn.dt[!(pas_alloc == "none" & act_alloc == "none"), unique(scenario)]

#' comes key'd
doses.dt <- readRDS(.args[2])[
  eval(datefilter)
][scenario %in% scns, .(
  scenario, realization, outcome, date, c.value
)]

intscns <- doses.dt[, unique(scenario)]

plt.dt <- setkeyv(
  doses.dt[scn.dt[scenario %in% intscns], on=.(scenario)],
  union(key(doses.dt), colnames(scn.dt))
)

rm(doses.dt)
gc()

plt.qs <- quantile(
  plt.dt,
  j = .(c.value), sampleby = "realization",
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
    breaks = c("risk", "ring", "none"),
    labels = c(
      ring="Ring Vaccination",
      risk="Risk-Based Strategy",
      none="Standard Only"
    ),
    aesthetics = c("color", "fill")
  ) +
  scale_y_continuous(
    name = sprintf(
      "Cumulative Vaccination"
    )
  ) +
  scale_x_null() +
  scale_linetype_manual("Vaccinate...", labels = c(conditional="Condtional on\nCase History", unconditional = "Unconditionally"), values = c(conditional="dashed", unconditional="solid")) +
  scale_alpha(range = c(0.01, 1)) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1), strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")
